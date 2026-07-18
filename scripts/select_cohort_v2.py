import sys, csv, json, glob, re
sys.path.insert(0,"src")
from vcf2report.phenopacket import load_phenopacket
SP="/private/tmp/claude-501/-Users-gbbarra-Documents-BbyB2/9c3eaae6-eebc-493f-a0b5-10b4b68c6a27/scratchpad"
def norm(c): return str(c).replace("chr","")

used_genes={r["gene"].upper() for r in csv.DictReader(open("data/synthetic_cohort/cohort.tsv"),delimiter="\t")}
used_coords={(norm(r["chrom"]),str(r["pos"]),r["ref"].upper(),r["alt"].upper()) for r in csv.DictReader(open("data/synthetic_cohort/cohort.tsv"),delimiter="\t")}
bg=[b for b in open(f"{SP}/bg_available.txt").read().split() if b]

def classify(hp, ref, alt):
    hp=hp or ""
    if "fs" in hp or "frameshift" in hp.lower(): return "frameshift_variant"
    if re.search(r"(Ter|\*)\)?$", hp) or hp.endswith("*"): return "stop_gained"
    if "Met1" in hp or "M1" in hp: return "start_lost"
    if re.search(r"del|dup|ins", hp):  # in-frame if len diff %3==0
        return "inframe_indel" if abs(len(ref)-len(alt))%3==0 and len(ref)!=len(alt) else "frameshift_variant"
    if re.search(r"p\.[A-Z][a-z][a-z]\d+[A-Z][a-z][a-z]", hp): return "missense_variant"
    if len(ref)==1 and len(alt)==1: return "missense_variant"  # SNV c.-level, assume missense-ish
    return "other"

# indexar phenopackets: 1 caso -> gene, variantes, hpo, disease
cases={}
for fp in glob.glob(f"{SP}/pp/**/*.json", recursive=True):
    try: d=load_phenopacket(fp)
    except Exception: continue
    hpo=d.get("hpo_terms",[])
    if len(hpo)<3: continue
    vs=[v for v in d.get("variants",[]) if all(v.get(k) for k in ("chrom","pos","ref","alt"))]
    if not vs: continue
    # variante primaria = a primeira; genotipo = todas
    for pv in vs:
        g=(pv.get("gene") or "").upper()
        if not g or g in used_genes: continue
        key=(norm(pv["chrom"]),str(pv["pos"]),pv["ref"].upper(),pv["alt"].upper())
        if key in used_coords: continue
        cons=classify(pv.get("hgvs_p"), pv["ref"], pv["alt"])
        # disease
        dis=d.get("disease") or (d.get("diseases") or [{}])
        disname=""
        cases.setdefault(g, {"gene":g,"chrom":"chr"+norm(pv["chrom"]),"pos":pv["pos"],"ref":pv["ref"].upper(),"alt":pv["alt"].upper(),
                             "cons":cons,"hpo":hpo,"vs":vs,"pv":pv})
        break

print(f"genes candidatos NOVOS (nao-usados, >=3 HPO): {len(cases)}")
from collections import Counter
print("distribuicao de consequencia disponivel:")
for c,n in Counter(v["cons"] for v in cases.values()).most_common(): print(f"  {n:4}  {c}")

# selecao com OVERSAMPLING de missense/inframe (VUS-producing)
target={"missense_variant":50,"inframe_indel":15,"stop_gained":13,"frameshift_variant":15,"start_lost":7}
bygroup={}
for v in cases.values(): bygroup.setdefault(v["cons"],[]).append(v)
sel=[]
for cons,want in target.items():
    pool=sorted(bygroup.get(cons,[]), key=lambda v:v["gene"])
    sel += pool[:want]
# completar 100 se faltar
if len(sel)<100:
    extra=[v for v in cases.values() if v not in sel]
    sel += extra[:100-len(sel)]
sel=sel[:100]
print(f"\nselecionados: {len(sel)}")
print("distribuicao selecionada:")
for c,n in Counter(v["cons"] for v in sel).most_common(): print(f"  {n:3}  {c}")
print(f"missense+inframe (VUS-producing): {sum(1 for v in sel if v['cons'] in ('missense_variant','inframe_indel'))}/100 (v1 era 45)")

# === atribuir backgrounds novos + emitir cohort_v2.tsv + plano fiel ===
def norm2(c): return str(c).replace("chr","")
rows=[]; plan={}
for i,v in enumerate(sel):
    sid=f"SYN-{101+i:03d}"; sample=bg[i]
    hpo=",".join(v["hpo"])
    rows.append([sid,sample,v["gene"],v["chrom"],str(v["pos"]),v["ref"],v["alt"],v["cons"],v.get("disease",""),hpo])
    # plano fiel: genotipo do phenopacket
    vs=v["vs"]; pv=v["pv"]
    pkey=(norm2(pv["chrom"]),str(pv["pos"]),pv["ref"].upper(),pv["alt"].upper())
    pz=next((x.get("zygosity") for x in vs if (norm2(x["chrom"]),str(x["pos"]),x["ref"].upper(),x["alt"].upper())==pkey),"het")
    others=[x for x in vs if (norm2(x["chrom"]),str(x["pos"]),x["ref"].upper(),x["alt"].upper())!=pkey]
    if (pz or "").startswith("hom"): plan[sid]={"mode":"hom"}
    elif others:
        o=others[0]; plan[sid]={"mode":"compound_het","chrom":"chr"+norm2(o["chrom"]),"pos2":int(o["pos"]),"ref2":o["ref"].upper(),"alt2":o["alt"].upper(),"zyg2":o.get("zygosity","het")}
    else: plan[sid]={"mode":"single_het"}

import csv as _csv
with open(f"{SP}/cohort_v2_101_200.tsv","w") as f:
    w=_csv.writer(f,delimiter="\t"); w.writerow(["syn_id","sample","gene","chrom","pos","ref","alt","consequence","disease","hpo"])
    for r in rows: w.writerow(r)
json.dump(plan, open(f"{SP}/v2_faithful_plan_101_200.json","w"), indent=0)
from collections import Counter as C
print(f"\n>>> cohort_v2_101_200.tsv escrito ({len(rows)} casos)")
print(f"    genotipos fieis: {dict(C(p['mode'] for p in plan.values()))}")
print(f"    backgrounds distintos dos usados: {len(set(r[1] for r in rows) - set(open(f'{SP}/used_bg.txt').read().split()))}/100")
print("    primeiros 3:")
for r in rows[:3]: print(f"      {r[0]} {r[1]} {r[2]} {r[3]}:{r[4]} {r[7]}")
