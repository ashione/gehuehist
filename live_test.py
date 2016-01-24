import os,sys
from scipy.stats import spearmanr as spr
import numpy as np
#SCIsdatabase = '~/IQA_CNN/SCIs/DistortedImages'
os.system('libsvm/svm-scale -l -1 -u 1 -s allrange live_train.txt > train_scale')
os.system('libsvm/svm-train  -s 3 -g 0.03 -c 1024 -b 1 -q train_scale allmodel')

if len(sys.argv)>1 and  sys.argv[1]=='retest':
    with open('live_test.txt') as fp :
        lines = fp.readlines()
        lines = map(lambda x : x.replace('\n',''),lines)
        image_info = map(lambda x : x.split(' '),lines)
    os.system('rm live_test_score.txt')
    cmds = map(lambda x : './brisquequality -im '+x[0]+' >> live_test_score.txt',image_info)
    for i,cm in enumerate(cmds):
        print i
        os.system(cm)
    #map(lambda x :os.system(x),cmds)

#image_name = map(lambda x : '_'.join(['cim'+str(x/49+1),str((x-x/49*49)/7+1),str((x-x/49*49)%7+1)])+'.bmp',range(0,980))
test_score = np.loadtxt('live_test_score.txt')
mos_score = np.loadtxt('live_mos.txt')
print spr(test_score,mos_score)
