#!/usr/bin/python3
import re
import os

#parser = argparse.ArgumentParser()
#parser.add_argument('-o', '--output', type=str, default='output.txt', help="Output file path")
#args = parser.parse_args()

#output_file = args.output

def blastfmt(f,f1):
#  output_filename = "raw_vdj.txt"
#  if os.path.exists(output_filename):
#    os.remove(output_filename)
  with open(f1, 'w') as g:
#    if os.path.exists(raw_vdj.txt):
#       os.remove(raw_vdj.txt)
    out0=('Query'+'\t'+'Gene'+'\t'+'Length'+'\t'+'Score'+'\t'+'Expect'+'\t'+'Identities'+'\t'+'Iden_percentage'+'\t'+'Gaps'+'\t'+'Gap_percentage'+'\t'+'Strand'+'\t'+'querystart'+'\t'+'queryend'+'\t'+'queryseq'+'\t'+'subjectstart'+'\t'+'subjectend'+'\t'+'subjectseq'+'\n')
    g.write(out0)
    sname = qseq = sseq = ""
    qstart = qend = 0
    sstart = send = 0
    le = 0
    score = expect = 0
    iden = idenper = gap = gapper = 0
    strand = 0
    #query='m64039_210630_173216/83953370/ccs'
    query=0
    for line in f:
         match = re.search("Query= .*", line)
         if match:
               q = (match.group())
               query = q.split(' ')[1]
             #  g.write(query)
         match = re.search(r">(.+)\n", line)
         if match:
               if qstart != 0:
                   out1=(str(query)+'\t'+sname+'\t'+le+'\t'+score+'\t'+expect+'\t'+iden+'\t'+idenper+'\t'+gap+'\t'+gapper+'\t'+strand+'\t'+qstart+'\t'+qend+'\t'+qseq+'\t'+sstart+'\t'+send+'\t'+sseq)
                   g.write(out1+'\n')
                   sname = qseq = sseq = ""
                   le = 0
                   score = expect = 0
                   iden = idenper = gap = gapper = 0
                   strand = 0
                   qstart = qend = 0
                   sstart = send = 0
              #     query=query1
                   sname = match.group(1)
               else:
                   sname = match.group(1)
               continue
         match = re.search(r"Query\s+(\d+)\s+(\S+)\s+(\d+)", line)
         if match:
               if qstart == 0:
                   qstart = match.group(1)
               qseq += match.group(2)
               qend = match.group(3)
               continue
         match = re.search(r"Sbjct\s+(\d+)\s+(\S+)\s+(\d+)", line)
         if match:
               if sstart == 0:
                  sstart = match.group(1)
               sseq += match.group(2)
               send = match.group(3)
               continue
         match = re.search("Length=\d+", line) 
         if match:
               length=(match.group()) 
               le = length.split('=')[1]
               continue
         match = re.search("Score = (\d+).*Expect = (\S+)", line)
         if match:
 #              if sstart == 0:
                   score = match.group(1)
                   expect = match.group(2)
         match = re.search("Identities = (\d+).* Gaps = (\S+).*",line)
         if match:
 #              if sstart == 0:
                   d = match.group()
                   dx = d.split('=',2)
                   dy = dx[1]
                   dy1 = dy.split(',',1)
                   d1 = dy1[0]
                   d2 = dx[2]
                   iden=d1.split('(',2)[0]
                   iden1=d1.split('(',2)[1]
                   idenper=iden1.split(')',2)[0]
                   gap=d2.split('(',2)[0]
                   gap1=d2.split('(',2)[1]
                   gapper=gap1.split(')',2)[0]
         match = re.search("Strand=\S+", line)
         if match:
               s = match.group()
               strand = s.split('=')[1]
         match = re.search("Effective search space used", line)
         if match:
              if qstart != 0:
                    out2=(str(query)+'\t'+sname+'\t'+le+'\t'+score+'\t'+expect+'\t'+iden+'\t'+idenper+'\t'+gap+'\t'+gapper+'\t'+strand+'\t'+qstart+"\t"+qend+"\t"+qseq+"\t"+sstart+"\t"+send+"\t"+sseq)
                    g.write(out2+'\n')
                    sname = qseq = sseq = ""
                    le = 0
                    score = expect = 0
                    iden = idenper = gap = gapper = 0
                    strand = 0
                    qstart = qend = 0
                    sstart = send = 0
    if qstart != 0:
          out3=(query+'\t'+sname+'\t'+le+'\t'+score+'\t'+expect+'\t'+iden+'\t'+idenper+'\t'+gap+'\t'+gapper+'\t'+strand+'\t'+qstart+"\t"+qend+"\t"+qseq+"\t"+sstart+"\t"+send+"\t"+sseq)
          g.write(out3+'\n')
    g.close()


#    g.seek(0)
#    lines = g.read()
        # return(line.strip())
#    return lines
    #with open("output.txt", 'r') as g:
     #     return g.read()
