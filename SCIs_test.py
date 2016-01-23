import os,sys
SCIsdatabase = '~/IQA_CNN/SCIs/DistortedImages'
image_name = map(lambda x : '_'.join(['cim'+str(x/49+1),str((x-x/49*49)/7+1),str((x-x/49*49)%7+1)])+'.bmp',range(0,980))
cmds = map(lambda x : './brisquequality -im '+os.path.join(SCIsdatabase,x),image_name)
map(lambda x :os.system(x),cmds)
