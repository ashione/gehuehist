import os,sys
from scipy.stats import spearmanr as spr
import numpy as np
#SCIsdatabase = '~/IQA_CNN/SCIs/DistortedImages'
for x in range(10,100):
    gamma = float(x)/1000
    os.system('libsvm/svm-scale -l -1 -u 1 -s allrange live_train.txt > train_scale')
    os.system('libsvm/svm-train  -s 3 -g {0} -c 2048 -b 1 -q train_scale allmodel'.format(gamma))
    #print str(x)

    if len(sys.argv)>1 and  sys.argv[1]=='retest':
        with open('live_test.txt') as fp :
            lines = fp.readlines()
            lines = map(lambda x : x.replace('\n',''),lines)
            image_info = map(lambda x : x.split(' '),lines)
            fp.close()
        os.system('rm live_test_score.txt')
        cmds = map(lambda x : './brisquequality -im '+x[0]+' >> live_test_score.txt',image_info)
        for i,cm in enumerate(cmds):
            print i
            os.system(cm)
        #map(lambda x :os.system(x),cmds)

    #image_name = map(lambda x : '_'.join(['cim'+str(x/49+1),str((x-x/49*49)/7+1),str((x-x/49*49)%7+1)])+'.bmp',range(0,980))
    test_score = np.loadtxt('live_test_score.txt')
    mos_score = np.loadtxt('live_mos.txt')
    sp_v  = spr(test_score,mos_score)
    print sp_v
    fp = open('live_svm_parms.txt','a')
    fp.write('{0} : spearmanr : {1}\n'.format(gamma,sp_v))
    fp.close()
