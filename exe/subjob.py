import sys, os
sys.path.append ('py')
import Submit as sb

os.system('rm para/*')
os.system('./gensublatt.exe '+sb.BASEPARA)
sb.makesubpara(sb.BASEPARA)
sb.submit()
