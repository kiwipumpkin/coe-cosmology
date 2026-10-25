cd /workspace/coepdf; timeout 150 tectonic main_rev5.tex > chk/build.log 2>&1; grep -i "^error" chk/build.log; N=$(pdfinfo main_rev5.pdf|awk '/Pages/{print $2}'); echo pages $N
for p in $(seq 1 $N); do echo "$p: $(pdftotext -f $p -l $p main_rev5.pdf - | grep -o "\(Figure\|Table\) [0-9]*:\|^Figure$" | tr "\n" " ") w=$(pdftotext -f $p -l $p main_rev5.pdf - | wc -w)"; done
