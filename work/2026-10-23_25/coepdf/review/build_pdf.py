import re,html
from reportlab.lib.pagesizes import A4,landscape
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import glob
f=glob.glob('/usr/share/fonts/**/DejaVuSans.ttf',recursive=True)[0]
pdfmetrics.registerFont(TTFont('DV',f));pdfmetrics.registerFont(TTFont('DVB',f.replace('DejaVuSans.ttf','DejaVuSans-Bold.ttf')))
from reportlab.pdfbase.pdfmetrics import registerFontFamily
registerFontFamily('DV',normal='DV',bold='DVB',italic='DV',boldItalic='DVB')
ss=getSampleStyleSheet()
B=ParagraphStyle('b',parent=ss['BodyText'],fontName='DV',fontSize=8.5,leading=11)
C=ParagraphStyle('c',parent=B,fontSize=7,leading=8.5)
H={1:ParagraphStyle('h1',parent=B,fontName='DVB',fontSize=14,leading=18,spaceBefore=6),2:ParagraphStyle('h2',parent=B,fontName='DVB',fontSize=11,leading=14,spaceBefore=6)}
def inl(t):
    t=html.escape(t,quote=False)
    t=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',t);t=re.sub(r'`(.+?)`',r'\1',t)
    return t
W=landscape(A4)[0]-60
def md(path,story):
    lines=open(path).read().split('\n');i=0
    while i<len(lines):
        l=lines[i]
        if l.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].startswith('|'):
                c=[x.strip().replace('\x00','|') for x in lines[i].replace('\\|','\x00').strip().strip('|').split('|')]
                if not all(re.fullmatch(r'-+',x) for x in c): rows.append(c)
                i+=1
            n=len(rows[0]);data=[[Paragraph(inl(x),C) for x in r[:n]]+[Paragraph('',C)]*(n-len(r)) for r in rows]
            lens=[min(max(max(len(r[j]) if j<len(r) else 0 for r in rows),12),400)**0.75 for j in range(n)];s=sum(lens)
            t=Table(data,colWidths=[W*x/s for x in lens],repeatRows=1)
            t.setStyle(TableStyle([('GRID',(0,0),(-1,-1),0.3,colors.grey),('BACKGROUND',(0,0),(-1,0),colors.lightgrey),('VALIGN',(0,0),(-1,-1),'TOP')]))
            story+= [t,Spacer(1,4)];continue
        m=re.match(r'(#+)\s+(.*)',l)
        if m: story.append(Paragraph(inl(m.group(2)),H[min(len(m.group(1)),2)]))
        elif l.strip(): story.append(Paragraph(('• '+inl(l.strip()[2:])) if l.strip().startswith('- ') else inl(l),B))
        i+=1
story=[Paragraph('Peer Review — Off-centre CoE cosmology (main_rev5)',H[1]),Paragraph('ARS academic-paper-reviewer workflow; NOT_CALIBRATED; single-executor sequential passes. Contents: editorial decision & top-10 issues; math-check table; five reports; Phase-0 field analysis.',B),Spacer(1,6)]
for p in ['07_editorial_decision.md','06_math_check.md','01_eic_journal_fit.md','02_methodology.md','03_domain.md','04_perspective.md','05_devils_advocate.md','00_field_analysis.md']:
    md(p,story);story.append(PageBreak())
SimpleDocTemplate('coe_review.pdf',pagesize=landscape(A4),leftMargin=30,rightMargin=30,topMargin=30,bottomMargin=30).build(story)
