import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class PROMPI_BURN():

    def __init__(self,filename):

        ff = open(filename,'r')
		
#        self.nt = int(ff.readline())
#        self.nnuc = int(ff.readline())

        header = ff.readline().split()
        self.nt = int(header[0])
        self.nnuc = int(header[1])
        self.tt = float(header[2])
        self.dd = float(header[3])
        print(self.nt,self.nnuc,self.tt,self.dd)

        ff.close()

        self.df = pd.read_csv(filename, skiprows=1,sep='\s+')
        
    def data(self,key):
	    return self.df[key]


class COCOCUBED():

    def __init__(self,filename1,filename2):

        self.dfe = pd.read_csv(filename1, skiprows=2,sep='\s+')
        self.dfx = pd.read_csv(filename2, skiprows=2,sep='\s+')
        self.filename1 = filename1
        self.filename2 = filename2

    def datae(self,key):
	    return self.dfe[key]        

    def datax(self,key):
	    return self.dfx[key]         


class COCOCUBED_TORCH():

    def __init__(self,filename1,filename2,filename3,filename4,filename5):

        self.dfe = pd.read_csv(filename1, skiprows=2,sep='\s+')
        self.dfx2 = pd.read_csv(filename2, skiprows=2,sep='\s+')
        self.dfx3 = pd.read_csv(filename3, skiprows=2,sep='\s+')
        self.dfx4 = pd.read_csv(filename4, skiprows=2,sep='\s+')
        self.dfx5 = pd.read_csv(filename5, skiprows=2,sep='\s+')
        
        self.filename1 = filename1
        self.filename2 = filename2

    def datae(self,key):
	    return self.dfe[key]        

    def datax2(self,key):
	    return self.dfx2[key]         

    def datax3(self,key):
	    return self.dfx3[key]

    def datax4(self,key):
	    return self.dfx4[key]         

    def datax5(self,key):
	    return self.dfx5[key]
        
class COCOCUBED_OBURN25():

    def __init__(self,filename1,filename2,filename3,filename4):

        self.dfe = pd.read_csv(filename1, skiprows=2,sep='\s+')
        self.dfx2 = pd.read_csv(filename2, skiprows=2,sep='\s+')
        self.dfx3 = pd.read_csv(filename3, skiprows=2,sep='\s+')
        self.dfx4 = pd.read_csv(filename4, skiprows=2,sep='\s+')
        
        self.filename1 = filename1
        self.filename2 = filename2

    def datae(self,key):
	    return self.dfe[key]        

    def datax2(self,key):
	    return self.dfx2[key]         

    def datax3(self,key):
	    return self.dfx3[key]

    def datax4(self,key):
	    return self.dfx4[key]         
