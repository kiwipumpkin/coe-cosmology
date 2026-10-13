import json
from jn_report import load,LAB
from views import sheet
SEL={"m":[("free",24),("cont",12),("loose",24),("web",24),("strict",24)],"n":[("free",32),("cont",20),("loose",32),("web",32),("strict",32)]}
for s,lst in SEL.items():
    items=[]
    for lv,N in lst:
        v=load(lv,s)[N]; r=v["show"]; tag="all rules met" if v["met"] else "rules not all met, best compromise"
        items.append((f"{LAB[lv]}\nN={N}, {tag}",r))
    sheet(items,("Mirror layouts allowed" if s=="m" else "No mirror layouts")+": best joint layouts at each rule level",f"side_views_{'mirror' if s=='m' else 'nomirror'}.png")
