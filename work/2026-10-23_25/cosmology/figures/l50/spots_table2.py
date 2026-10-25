import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
R=[('A', '1', 'CMB dipole', '84°, −48°', 'cold', '< 0', '+2.8%', '✗', 'Our motion → distance/velocity'), ('A', '2', 'Bulk flow (W09)', '107°, −8°', 'cold', '< 0', '−1.6%', '✓', 'Recession velocity → distance'), ('A', '3', 'CMB Cold Spot', '209°, −57°', 'cold', '< 0', '−2.9%', '✓', 'Eridanus supervoid → distances'), ('A', '4', 'CMB dipole', '264°, +48°', 'hot', '> 0', '+0.2%', '✓', 'Our motion → distance/velocity'), ('A', '5', 'Bulk flow (W09)', '287°, +8°', 'hot', '> 0', '+2.5%', '✓', 'Recession velocity → distance'), ('B', '1', 'NVSS dipole', '16°, −46°', 'cold', '< 0', 'low zone (crossing 14°)', '✓', 'Source density → sample composition'), ('B', '2', 'Hemispherical asym.', '50°, +20°', 'cold', '< 0', '−4.4%', '✓', 'Power direction → sample/foreground'), ('B', '3', 'CatWISE dipole', '58°, −31°', 'cold', '< 0', '−4.9%', '✓', 'Matter dipole → sample composition'), ('B', '4', 'NVSS dipole', '196°, +46°', 'hot', '> 0', '+1.3%', '✓', 'Source density → sample composition'), ('B', '5', 'Hemispherical asym.', '230°, −20°', 'hot', '> 0', '+2.8%', '✓', 'Power direction → sample/foreground'), ('B', '6', 'CatWISE dipole', '238°, +31°', 'hot', '> 0', '+2.1%', '✓', 'Matter dipole → sample composition'), ('B', '7', 'Loop I / NPS', '329°, +17.5°', 'hot', '> 0', '+0.9%', '✓', 'Soft X-ray foreground + N_H → tilts B')]
H=['Tied to','#','Spot','ℓ, b','Hot/cold','Needs','Zone / residual (±25°)','Fits?','Why A or B']
fig,ax=plt.subplots(figsize=(17,5));ax.axis('off')
t=ax.table(cellText=R,colLabels=H,loc='center',cellLoc='center',colWidths=[.06,.03,.16,.09,.08,.07,.05,.05,.41])
t.auto_set_font_size(False);t.auto_set_column_width(list(range(9)));t.set_fontsize(10);t.scale(1,1.6)
for (r,c),cl in t.get_celld().items():
  if r==0: cl.set_facecolor('#e8e8e8');cl.get_text().set_weight('bold')
  elif c==4: v=R[r-1][4];cl.get_text().set_color({'hot':'red','cold':'blue'}.get(v,'0.4'));cl.get_text().set_weight('bold' if v in('hot','cold') else 'normal')
  elif c==7 and R[r-1][7] in '✓✗': cl.get_text().set_color('green' if R[r-1][7]=='✓' else 'red')
ax.set_title('Hot spot needs residual > 0, cold spot needs residual < 0 (Migkas-fitted regions excluded)\nA = normalization ↔ distance/velocity (H₀ ∝ A^½);   B = slope ↔ cluster physics, foreground, sample',fontsize=11)
fig.text(.5,.02,'Score: A 4/5, B 7/7, total 11/12 matches',ha='center',weight='bold',fontsize=12)
fig.savefig('spots_table2.png',dpi=150,bbox_inches='tight',facecolor='white')
