#Quintin Tyrell Davis
#Comments 2/5/2013
#This code was used as the basis for a presentation at the 2011 Wind River Conference On Prokaryotic Biology, June 2011, Estes Park Colorado. 
#


class matrix:
    def __init__(self,cont = [[]]):
        self.m = cont
        
    def printMat(self):
        i = 0
        while (i < len(self.m)):
            print(self.m[i])
            
    def add(self, other):
        i = 0
        k = 0
        result = matrix([[]])
        #check matrix dimensions for agreement
        if (len(self.m) != len(other.m)):
            print("Error: matrix dimensions do not agree")
            return 0
        while (i<len(self.m)):
            if (len(self.m[i]) != len(other.m[i])):
                print("Error: matrix dimensions do not match")
                return 0    
            i += 1
        #add the matrices
        i = 0
        while (i < len(self.m)):
            k = 0
            while (k < len(self.m[i])):
                if ((i+1) == len(result.m)):
                    #print(result.m, result.m[i])
                    result.m[i].append( self.m[i][k] + other.m[i][k] )
                    
                else:
                    result.m.append([self.m[i][k] + other.m[i][k]])
                k += 1
            i += 1    
        #result.display()      #debug statement to print result       
        return result
                
        
    def times(self, other):
        if len(self.m[0]) != len(other.m):
            print("Error: n & m dimensions do not match")
            return 0
        result = [[]]
        temp = 0
        i = 0
        while (i < len(self.m)):
            c = 0
            while (c < len(other.m[0])):
                   k = 0
                   temp = 0
                   while (k < len(other.m)):
                       temp += (self.m[i][k]) * (other.m[k][c])
                       k += 1 
                   if (len(result) == (i + 1)):
                       result[i] += [temp]
                   else:
                       result.append([temp])
                   c += 1
            i += 1
                
        return matrix(result)
    
    def sTimes(self, scalar):
        result = []#matrix(self.m)
        i = 0
        while (i < len(self.m)):
            k = 0
            temp = []
            while (k < len(self.m[0])):
                temp.append((self.m[i][k]) * (scalar))
                k += 1
            result.append(temp)    
            i += 1
           
            
            
        print (result)
        return matrix(result)

    def display(self):
        i = 0
        k = 0
        while (i < len(self.m)):
            k = 0
            while (k < len(self.m[i])):
                print(repr(self.m[i][k]).rjust(4), end = '  ')
                k += 1
            print('')
            i += 1

    def transp(self):
        """Transpose matrix"""


class network:
    
    def __init__(self, fIn = 0, runCom = [0,1.0,0.01] ):
        if (fIn == 0):
            print("Error, no input file read""")
            return 0
        self.vEnz = fIn[0]
        self.vMet = fIn[1]
        self.mStoch = fIn[2]
        self.mKin = fIn[3]
        self.mAdj = []
        self.vENames = []
        self.vMetNames = []
        self.runCom = runCom
 
    def tRun(self, file):#, tDur = 1.0, tRes = 0.01):
        """variable timestep kinetics simulation"""
        #file = openResFile()
        tDur = self.runCom[1]
        tRes = self.runCom[2]
        mtMetResults = self.vMet
        mtEnzResults = self.vEnz
        i = 0
        k = 0
        writeResFile(file,str(i))
        while (k < len(self.vMet.m)):
            writeResFile(file,str('\t') +  str(self.vMet.m[k][0]))
            k+=1
        writeResFile(file,'\n')
                     
        while ( i < tDur ):
            mtMetResults = self.mStoch.times(self.buildFluxVector().sTimes(tRes))
            x = 0
            #self.vMet = self.vMet.add(mtMetResults)
            while ( x < len(self.vMet.m)):
                self.vMet.m[x][0] = max((self.vMet.m[x][0] + mtMetResults.m[x][0]),0.0)            #Right now just modify state of vMet
                x += 1
            i += tRes
            k = 0
            writeResFile(file,str(i))
            while (k < len(self.vMet.m)):
                writeResFile(file,str('\t') +  str(self.vMet.m[k][0]))         #0.0 max statement is a fudge factor for now to avoid problems arising from too large of step resolution
                k+=1
            writeResFile(file,'\n')

    def paraEnzRun(self, file, Enz = 0, rng = [0.0,1.0], rngRes = 0.1):
        cnc = rng[0]
        f = readFile()
        metVec = getMetVec(f)

        while (cnc < rng[1]):
            self.vMet = matrix(getConcs(metVec))
            A.vEnz.m[Enz][0] = cnc
            writeResFile(file,'\n')#conc\t'+str(cnc)+'\n')
            self.tRun(file)
            cnc += rngRes

    def phasePlane(self, file, Met = 0, rng = [0.0,1.0], rngRes = 0.1, consts = [], const = []):
               
        cnc = rng[0]
        while (cnc < rng[1]):
            self.vMet = self.vMet.sTimes(0)
            i = 0
            while (i < len(consts)):
                self.vMet.m[consts[i]][0] = const[i]
                i += 1
            self.vMet.m[Met][0] = cnc
            writeResFile(file,'\n')#+'conc\t' + str(cnc) + '\n')
            self.tRun(file)
            cnc += rngRes
            
    def paraMetRun(self, file, Met = 0, rng = [0.0,1.0], rngRes = 0.1):
        f = readFile()
        metVec = getMetVec(f)
        #netList = getNetList(f)
        #runCom = getRunCom(f)
        #self.vEnz = matrix(getConcs(netList))
        #self.vMet = matrix(getConcs(metVec))
        #self.mKConsts = matrix(buildKinConstMat(netList))
        #print(vEnz,vMet,mKConsts)
        #nodes = getNodeArray(netList)
        #print(nodes)
        #mStoch = getStochMat(nodes)
        #mStoch.display()
        #runCom = parseRunCom(getRunCom(f))
        
        cnc = rng[0]
        while (cnc < rng[1]):
            self.vMet = matrix(getConcs(metVec))
            self.vMet.m[Met][0] = cnc
            writeResFile(file,'\n')#+'conc\t' + str(cnc) + '\n')
            self.tRun(file)
            cnc += rngRes

    def detRxnTypes(self,i = 0):
        s = 0
        p = 0
        
        #"""Below calculates the number of reactants by counting down the columns of the stochiometric matrix"""
        #while (i < len(self.mStoch.m[0])):
        k = 0
        while (k < len(self.mStoch.m)):
            """if (self.mStoch.m[k][i] == 0):
                s = s #do nothing"""
            if (self.mStoch.m[k][i] < 0):
                s += 1
            if (self.mStoch.m[k][i] > 0):
                p += 1
            k += 1
            #i += 1
        if ((s == 1) and (p == 1)):
            #print('rxn is uni-uni')
            return 0
        elif ((s == 1) and (p == 2)):
            #print('rxn is uni-bi')
            return 1
        elif ((s == 2) and (p == 1)):
            #print('rxn is bi-uni')
            return 2
        elif ((s == 2) and (p == 2)):
            #print('rxn is bi-bi')
            return 3
    
#    def buildFluxVecComp(self, Enz = 0, typ = 0):
#        fluxComp = 0.0
#        i = 0
#        n = []
#        out = []
#        while (i < len(self.vMet.m)):
#            if (self.mStoch.m[i][Enz] < 0):
#                n.append(i)
#            elif (self.mStoch.m[i][Enz] > 0):
#                out.append(i)
#            i += 1
#        #print('in = ', n, 'out = ', out, 'enzyme = ' , Enz)#debug statement
#        if ( typ == 0 ):
#            #The equation for uni-uni reactant kinetics is from Enzyme Kinetics and Mechanism, Cook and Cleland 2007 p 15
#            fluxComp = ((self.vEnz.m[Enz][0] * ((self.mKin.m[Enz][1] * self.vMet.m[n[0]][0])-(self.mKin.m[Enz][3]*((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*self.vMet.m[out[0]][0])))/((self.mKin.m[Enz][0])+(self.vMet.m[n[0]][0])+(((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*self.vMet.m[out[0]][0])))
#            #elif (i > self.mStoch.m[i][Enz] > 0):
#        if ( typ == 1 ):
#            #the equation for uni-bi reactant kinetics is from Enzyme Kinetics and Mechanism, Cook and Cleland 2007 p 15
#            fluxComp = ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * self.vMet.m[n[0]][0])-   ( self.mKin.m[Enz][3] * self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + self.vMet.m[n[0]][0] + ((self.mKin.m[Enz][0] * self.vMet.m[out[0]][0])/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * self.vMet.m[out[1]][0])/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*self.vMet.m[out[0]][0]*self.vMet.m[out[1]][0])/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
#        if ( typ == 2 ):
#            #for bi-uni rxn, ibid
#            fluxComp = - ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * self.vMet.m[n[0]][0])-   ( self.mKin.m[Enz][3] * self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + self.vMet.m[n[0]][0] + ((self.mKin.m[Enz][0] * self.vMet.m[out[0]][0])/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * self.vMet.m[out[1]][0])/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*self.vMet.m[out[0]][0]*self.vMet.m[out[1]][0])/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
#        if ( typ == 3 ):
#            #bi-bi rxn, 
#            #fluxComp = ((self.vEnz.m[Enz][0] * (self.mKin.m[Enz][1] * ( ( self.vMet.m[n[0]][0] * self.vMet.m[n[1]][0] )  / ( self.mKin.m[Enz][0] * self.mKin.m[Enz][2] ) ) - self.mKin.m[Enz][3] * ( (self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] ) / ( self.mKin.m[Enz][4] * self.mKin.m[Enz][5] ) ) ) ) / ( ( 1 + (self.vMet.m[n[0]][0] / self.mKin.m[Enz][0] ) + (self.vMet.m[out[0]][0] / self.mKin.m[Enz][4] ) ) * ( 1 + ( self.vMet.m[n[1]][0] / self.mKin.m[Enz][2] ) + (self.vMet.m[out[1]][0] / self.mKin.m[Enz][5] ) ) ) )
#            #from 221
#            fluxComp = ( self.vEnz.m[Enz][0] *(self.mKin.m[Enz][1]*((self.vMet.m[n[0]][0] * self.vMet.m[n[1]][0] ) / ( self.mKin.m[Enz][0] * self.mKin.m[Enz][2] ) )-(self.mKin.m[Enz][3] * (self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] )/(self.mKin.m[Enz][4] * self.mKin.m[Enz][5])))) / ((1+(self.vMet.m[n[0]][0] / self.mKin.m[Enz][0] )+(self.vMet.m[out[0]][0] / self.mKin.m[Enz][4]))*(1+(self.vMet.m[n[1]][0] / self.mKin.m[Enz][2] )+( self.vMet.m[out[1]][0] / self.mKin.m[Enz][5])))
#
#        return fluxComp
#

    def buildFluxVecComp(self, Enz = 0, typ = 0):
        #Runge-Kutta 2/ improved Euler's
        fluxComp = 0.0
        fluxComp1 = 0.0
        fluxComp2 = 0.0
        i = 0
        n = []
        out = []
        while (i < len(self.vMet.m)):
            if (self.mStoch.m[i][Enz] < 0):
                n.append(i)
            elif (self.mStoch.m[i][Enz] > 0):
                out.append(i)
            i += 1
        #print('in = ', n, 'out = ', out, 'enzyme = ' , Enz)#debug statement
        if ( typ == 0 ):
            #"""The equation for uni-uni reactant kinetics is from Enzyme Kinetics and Mechanism, Cook and Cleland 2007 p 15""""
            fluxComp1 = ((self.vEnz.m[Enz][0] * ((self.mKin.m[Enz][1] * self.vMet.m[n[0]][0])-(self.mKin.m[Enz][3]*((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*self.vMet.m[out[0]][0])))/((self.mKin.m[Enz][0])+(self.vMet.m[n[0]][0])+(((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*self.vMet.m[out[0]][0])))
            fluxComp2 = ((self.vEnz.m[Enz][0] * ((self.mKin.m[Enz][1] * (self.vMet.m[n[0]][0] - (fluxComp1 * self.runCom[2])))-(self.mKin.m[Enz][3]*((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*(self.vMet.m[out[0]][0] + (fluxComp1 * self.runCom[2])))))/((self.mKin.m[Enz][0])+((self.vMet.m[n[0]][0] - (fluxComp1 * self.runCom[2])))+(((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*(self.vMet.m[out[0]][0] + (fluxComp1 * self.runCom[2])))))
            #elif (i > self.mStoch.m[i][Enz] > 0):
        if ( typ == 1 ):
            #"""The equation for uni-bi reactant kinetics is from Enzyme Kinetics and Mechanism, Cook and Cleland 2007 p 15""""
            fluxComp1 = ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * self.vMet.m[n[0]][0])-   ( self.mKin.m[Enz][3] * self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + self.vMet.m[n[0]][0] + ((self.mKin.m[Enz][0] * self.vMet.m[out[0]][0])/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * self.vMet.m[out[1]][0])/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*self.vMet.m[out[0]][0]*self.vMet.m[out[1]][0])/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
            fluxComp2 = ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * (self.vMet.m[n[0]][0] - (fluxComp1 * self.runCom[2])))-   ( self.mKin.m[Enz][3] * (self.vMet.m[out[0]][0] + (fluxComp1 * self.runCom[2])) * (self.vMet.m[out[1]][0] + (fluxComp1 * self.runCom[2])) * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + (self.vMet.m[n[0]][0] - (fluxComp1 * self.runCom[2])) + ((self.mKin.m[Enz][0] * (self.vMet.m[out[0]][0] + (fluxComp1 * self.runCom[2])))/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * (self.vMet.m[out[1]][0] + (fluxComp1 * self.runCom[2])))/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*(self.vMet.m[out[0]][0] + (fluxComp1 * self.runCom[2]))*(self.vMet.m[out[1]][0] + (fluxComp1 * self.runCom[2])))/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
        if ( typ == 2 ):
            #for bi-uni rxn, ibid
            fluxComp1 = - ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * self.vMet.m[n[0]][0])-   ( self.mKin.m[Enz][3] * self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + self.vMet.m[n[0]][0] + ((self.mKin.m[Enz][0] * self.vMet.m[out[0]][0])/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * self.vMet.m[out[1]][0])/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*self.vMet.m[out[0]][0]*self.vMet.m[out[1]][0])/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
            fluxComp2 = - ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * (self.vMet.m[n[0]][0] - (fluxComp1 * self.runCom[2])))-   ( self.mKin.m[Enz][3] * (self.vMet.m[out[0]][0] + (fluxComp1 * self.runCom[2])) * (self.vMet.m[out[1]][0] + (fluxComp1 * self.runCom[2])) * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + (self.vMet.m[n[0]][0] - (fluxComp1 * self.runCom[2])) + ((self.mKin.m[Enz][0] * (self.vMet.m[out[0]][0] + (fluxComp1 * self.runCom[2])))/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * (self.vMet.m[out[1]][0] + (fluxComp1 * self.runCom[2])))/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*(self.vMet.m[out[0]][0] + (fluxComp1 * self.runCom[2]))*(self.vMet.m[out[1]][0] + (fluxComp1 * self.runCom[2])))/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
        if ( typ == 3 ):
            #bi-bi rxn, 
            #fluxComp = ((self.vEnz.m[Enz][0] * (self.mKin.m[Enz][1] * ( ( self.vMet.m[n[0]][0] * self.vMet.m[n[1]][0] )  / ( self.mKin.m[Enz][0] * self.mKin.m[Enz][2] ) ) - self.mKin.m[Enz][3] * ( (self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] ) / ( self.mKin.m[Enz][4] * self.mKin.m[Enz][5] ) ) ) ) / ( ( 1 + (self.vMet.m[n[0]][0] / self.mKin.m[Enz][0] ) + (self.vMet.m[out[0]][0] / self.mKin.m[Enz][4] ) ) * ( 1 + ( self.vMet.m[n[1]][0] / self.mKin.m[Enz][2] ) + (self.vMet.m[out[1]][0] / self.mKin.m[Enz][5] ) ) ) )
            #from 221
            fluxComp1 = ( self.vEnz.m[Enz][0] *(self.mKin.m[Enz][1]*((self.vMet.m[n[0]][0] * self.vMet.m[n[1]][0] ) / ( self.mKin.m[Enz][0] * self.mKin.m[Enz][2] ) )-(self.mKin.m[Enz][3] * (self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] )/(self.mKin.m[Enz][4] * self.mKin.m[Enz][5])))) / ((1+(self.vMet.m[n[0]][0] / self.mKin.m[Enz][0] )+(self.vMet.m[out[0]][0] / self.mKin.m[Enz][4]))*(1+(self.vMet.m[n[1]][0] / self.mKin.m[Enz][2] )+( self.vMet.m[out[1]][0] / self.mKin.m[Enz][5])))
            fluxComp2 = ( self.vEnz.m[Enz][0] *(self.mKin.m[Enz][1]*(((self.vMet.m[n[0]][0] - (fluxComp1 * self.runCom[2])) * (self.vMet.m[n[1]][0] - (fluxComp1 * self.runCom[2])) ) / ( self.mKin.m[Enz][0] * self.mKin.m[Enz][2] ) )-(self.mKin.m[Enz][3] * ((self.vMet.m[out[0]][0] + (fluxComp1 * self.runCom[2])) * (self.vMet.m[out[1]][0] + (fluxComp1 * self.runCom[2])) )/(self.mKin.m[Enz][4] * self.mKin.m[Enz][5])))) / ((1+((self.vMet.m[n[0]][0] - (fluxComp1 * self.runCom[2])) / self.mKin.m[Enz][0] )+((self.vMet.m[out[0]][0] + (fluxComp1 * self.runCom[2])) / self.mKin.m[Enz][4]))*(1+((self.vMet.m[n[1]][0] - (fluxComp1 * self.runCom[2])) / self.mKin.m[Enz][2] )+( (self.vMet.m[out[1]][0] + (fluxComp1 * self.runCom[2])) / self.mKin.m[Enz][5])))
        fluxComp = (fluxComp1 + fluxComp2 ) / 2
        return fluxComp

    def buildFluxVecCompRK4(self, Enz = 0, typ = 0):
        #Runge-Kutta 4th order method (RK4)
        i = 0
        n = []
        out = []
        while (i < len(self.vMet.m)):
            if (self.mStoch.m[i][Enz] < 0):
                n.append(i)
            elif (self.mStoch.m[i][Enz] > 0):
                out.append(i)
            i += 1
        #print('in = ', n, 'out = ', out, 'enzyme = ' , Enz)#debug statement
        if ( typ == 0 ):
            #"""The equation for uni-uni reactant kinetics is from Enzyme Kinetics and Mechanism, Cook and Cleland 2007 p 15""""
            fluxComp1 = ((self.vEnz.m[Enz][0] * ((self.mKin.m[Enz][1] * self.vMet.m[n[0]][0])-(self.mKin.m[Enz][3]*((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*self.vMet.m[out[0]][0])))/((self.mKin.m[Enz][0])+(self.vMet.m[n[0]][0])+(((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*self.vMet.m[out[0]][0])))
            fluxComp2 = ((self.vEnz.m[Enz][0] * ((self.mKin.m[Enz][1] * ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp1 ) ))-(self.mKin.m[Enz][3]*((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp1 ) ))))/((self.mKin.m[Enz][0])+(( self.vMet.m[n[0]][0] - ( 0.5*fluxComp1 ) ))+(((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp1 ) ))))
            fluxComp3 = ((self.vEnz.m[Enz][0] * ((self.mKin.m[Enz][1] * ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp2 ) ))-(self.mKin.m[Enz][3]*((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp2 ) ))))/((self.mKin.m[Enz][0])+(( self.vMet.m[n[0]][0] - ( 0.5*fluxComp2 ) ))+(((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp2 ) ))))
            fluxComp4 = ((self.vEnz.m[Enz][0] * ((self.mKin.m[Enz][1] * ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp3 ) ))-(self.mKin.m[Enz][3]*((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp3 ) ))))/((self.mKin.m[Enz][0])+(( self.vMet.m[n[0]][0] - ( 0.5*fluxComp3 ) ))+(((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2]))*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp3 ) ))))
        if ( typ == 1 ):
            #"""The equation for uni-bi reactant kinetics is from Enzyme Kinetics and Mechanism, Cook and Cleland 2007 p 15""""
            fluxComp1 = ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * self.vMet.m[n[0]][0])-   ( self.mKin.m[Enz][3] * self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + self.vMet.m[n[0]][0] + ((self.mKin.m[Enz][0] * self.vMet.m[out[0]][0])/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * self.vMet.m[out[1]][0])/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*self.vMet.m[out[0]][0]*self.vMet.m[out[1]][0])/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
            fluxComp2 = ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp1 ) ))-   ( self.mKin.m[Enz][3] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp1 ) ) * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp1 ) ) * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp1 ) ) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp1 ) ))/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp1 ) ))/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp1 ) )*( self.vMet.m[out[1]][0] + ( 0.5*fluxComp1 ) ))/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
            fluxComp3 = ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp2 ) ))-   ( self.mKin.m[Enz][3] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp2 ) ) * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp2 ) ) * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp2 ) ) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp2 ) ))/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp2 ) ))/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp2 ) )*( self.vMet.m[out[1]][0] + ( 0.5*fluxComp2 ) ))/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
            fluxComp4 = ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp3 ) ))-   ( self.mKin.m[Enz][3] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp3 ) ) * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp3 ) ) * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp3 ) ) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp3 ) ))/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp3 ) ))/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp3 ) )*( self.vMet.m[out[1]][0] + ( 0.5*fluxComp3 ) ))/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
        if ( typ == 2 ):
           fluxComp1 = - ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * self.vMet.m[n[0]][0])-   ( self.mKin.m[Enz][3] * self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + self.vMet.m[n[0]][0] + ((self.mKin.m[Enz][0] * self.vMet.m[out[0]][0])/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * self.vMet.m[out[1]][0])/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*self.vMet.m[out[0]][0]*self.vMet.m[out[1]][0])/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
           fluxComp2 = - ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp1 ) ))-   ( self.mKin.m[Enz][3] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp1 ) ) * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp1 ) ) * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp1 ) ) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp1 ) ))/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp1 ) ))/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp1 ) )*( self.vMet.m[out[1]][0] + ( 0.5*fluxComp1 ) ))/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
           fluxComp3 = - ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp2 ) ))-   ( self.mKin.m[Enz][3] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp2 ) ) * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp2 ) ) * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp2 ) ) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp2 ) ))/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp2 ) ))/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp2 ) )*( self.vMet.m[out[1]][0] + ( 0.5*fluxComp2 ) ))/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
           fluxComp4 = - ( (self.vEnz.m[Enz][0] *((self.mKin.m[Enz][1] * ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp3 ) ))-   ( self.mKin.m[Enz][3] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp3 ) ) * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp3 ) ) * ((self.mKin.m[Enz][0])/(self.mKin.m[Enz][2] * self.mKin.m[Enz][4])) ))) / (self.mKin.m[Enz][0] + ( self.vMet.m[n[0]][0] - ( 0.5*fluxComp3 ) ) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[0]][0] + ( 0.5*fluxComp3 ) ))/(self.mKin.m[Enz][2])) + ((self.mKin.m[Enz][0] * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp3 ) ))/(self.mKin.m[Enz][4])) + ((self.mKin.m[Enz][0]*( self.vMet.m[out[0]][0] + ( 0.5*fluxComp3 ) )*( self.vMet.m[out[1]][0] + ( 0.5*fluxComp3 ) ))/(self.mKin.m[Enz][2]*self.mKin.m[Enz][4]))))
        if ( typ == 3 ):
            #bi-bi rxn, 
            #fluxComp = ((self.vEnz.m[Enz][0] * (self.mKin.m[Enz][1] * ( ( self.vMet.m[n[0]][0] * self.vMet.m[n[1]][0] )  / ( self.mKin.m[Enz][0] * self.mKin.m[Enz][2] ) ) - self.mKin.m[Enz][3] * ( (self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] ) / ( self.mKin.m[Enz][4] * self.mKin.m[Enz][5] ) ) ) ) / ( ( 1 + (self.vMet.m[n[0]][0] / self.mKin.m[Enz][0] ) + (self.vMet.m[out[0]][0] / self.mKin.m[Enz][4] ) ) * ( 1 + ( self.vMet.m[n[1]][0] / self.mKin.m[Enz][2] ) + (self.vMet.m[out[1]][0] / self.mKin.m[Enz][5] ) ) ) )
            #from 221
           fluxComp1 = ( self.vEnz.m[Enz][0] *(self.mKin.m[Enz][1]*((self.vMet.m[n[0]][0] * self.vMet.m[n[1]][0] ) / ( self.mKin.m[Enz][0] * self.mKin.m[Enz][2] ) )-(self.mKin.m[Enz][3] * (self.vMet.m[out[0]][0] * self.vMet.m[out[1]][0] )/(self.mKin.m[Enz][4] * self.mKin.m[Enz][5])))) / ((1+(self.vMet.m[n[0]][0] / self.mKin.m[Enz][0] )+(self.vMet.m[out[0]][0] / self.mKin.m[Enz][4]))*(1+(self.vMet.m[n[1]][0] / self.mKin.m[Enz][2] )+( self.vMet.m[out[1]][0] / self.mKin.m[Enz][5])))
           fluxComp2 = ( self.vEnz.m[Enz][0] *(self.mKin.m[Enz][1]*((( self.vMet.m[n[0]][0] - ( 0.5*fluxComp1 ) ) * ( self.vMet.m[n[1]][0] - ( 0.5*fluxComp1 ) ) ) / ( self.mKin.m[Enz][0] * self.mKin.m[Enz][2] ) )-(self.mKin.m[Enz][3] * (( self.vMet.m[out[0]][0] + ( 0.5*fluxComp1 ) ) * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp1 ) ) )/(self.mKin.m[Enz][4] * self.mKin.m[Enz][5])))) / ((1+(( self.vMet.m[n[0]][0] - ( 0.5*fluxComp1 ) ) / self.mKin.m[Enz][0] )+(( self.vMet.m[out[0]][0] + ( 0.5*fluxComp1 ) ) / self.mKin.m[Enz][4]))*(1+(( self.vMet.m[n[1]][0] - ( 0.5*fluxComp1 ) ) / self.mKin.m[Enz][2] )+( ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp1 ) ) / self.mKin.m[Enz][5])))
           fluxComp3 = ( self.vEnz.m[Enz][0] *(self.mKin.m[Enz][1]*((( self.vMet.m[n[0]][0] - ( 0.5*fluxComp2 ) ) * ( self.vMet.m[n[1]][0] - ( 0.5*fluxComp2 ) ) ) / ( self.mKin.m[Enz][0] * self.mKin.m[Enz][2] ) )-(self.mKin.m[Enz][3] * (( self.vMet.m[out[0]][0] + ( 0.5*fluxComp2 ) ) * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp2 ) ) )/(self.mKin.m[Enz][4] * self.mKin.m[Enz][5])))) / ((1+(( self.vMet.m[n[0]][0] - ( 0.5*fluxComp2 ) ) / self.mKin.m[Enz][0] )+(( self.vMet.m[out[0]][0] + ( 0.5*fluxComp2 ) ) / self.mKin.m[Enz][4]))*(1+(( self.vMet.m[n[1]][0] - ( 0.5*fluxComp2 ) ) / self.mKin.m[Enz][2] )+( ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp2 ) ) / self.mKin.m[Enz][5])))
           fluxComp4 = ( self.vEnz.m[Enz][0] *(self.mKin.m[Enz][1]*((( self.vMet.m[n[0]][0] - ( 0.5*fluxComp3 ) ) * ( self.vMet.m[n[1]][0] - ( 0.5*fluxComp3 ) ) ) / ( self.mKin.m[Enz][0] * self.mKin.m[Enz][2] ) )-(self.mKin.m[Enz][3] * (( self.vMet.m[out[0]][0] + ( 0.5*fluxComp3 ) ) * ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp3 ) ) )/(self.mKin.m[Enz][4] * self.mKin.m[Enz][5])))) / ((1+(( self.vMet.m[n[0]][0] - ( 0.5*fluxComp3 ) ) / self.mKin.m[Enz][0] )+(( self.vMet.m[out[0]][0] + ( 0.5*fluxComp3 ) ) / self.mKin.m[Enz][4]))*(1+(( self.vMet.m[n[1]][0] - ( 0.5*fluxComp3 ) ) / self.mKin.m[Enz][2] )+( ( self.vMet.m[out[1]][0] + ( 0.5*fluxComp3 ) ) / self.mKin.m[Enz][5])))
        fluxComp = (fluxComp1 + 2*fluxComp2 + 2* fluxComp3 + fluxComp4) / 6
        return fluxComp



    def buildFluxVector(self):
        vFlux = []
        #tempFlux = 0
        k = 0
        while(k < len(self.vEnz.m)):
            i = 0
            tempFlux = 0
            tempFlux = self.buildFluxVecComp(k,self.detRxnTypes(k))
            vFlux.append([tempFlux])
            k += 1
        return matrix(vFlux)

    def buildFluxVector0(self): #deprecated 5/19/201
        vFlux = []
        #tempFlux = 0
        k = 0
        while(k < len(self.vEnz.m)):
            i = 0
            tempFlux = 0
            
            while (i < len(self.vMet.m)):
                if (self.mStoch.m[i][k] < 0):
                    tempFlux += (self.vEnz.m[k][0] * self.mKin.m[k][1] * self.vMet.m[i][0]) / (self.mKin.m[k][0] + self.vMet.m[i][0])
                elif (i > self.mStoch.m[i][k] > 0):
                    tempFlux -= (self.vEnz.m[k][0] *  self.mKin.m[k][3] * self.vMet.m[i][0] ) / (self.mKin.m[k][2] + self.vMet.m[i][0])
                i +=1
                #print('i = ', i, 'k = ', k )
            k += 1
            vFlux.append([tempFlux])    
        #print(vFlux)
        return matrix(vFlux)


    def ppRun(self):
        """run kinetics simulation across a concentration sweep for two metabolites, output phase-plane plot data"""

    def paraRun(self):
        """run kinetics simulation with parameter sweep"""

    def openWFile(self):
        """Open results file for writing to"""
        
    def showStartConditions(self):
        print("\nEnzyme initial concentration = ")
        self.vEnz.display()
        print("\nMetabolite initial concentration = \n")
        self.vMet.display()
        print("\nStochiometric matrix = \n")
        self.mStoch.display()
        print("\nKinetic parameters = \n")
        self.mKin.display()
        
def openResFile():
    import time
    path = '/Volumes/QUINTIN/Project/MetabolSiM/Results/MSResults' + time.asctime(time.localtime()).replace(' ','').replace(':','') + '.txt'
    f = open(path, 'w')
    return f
    
def writeResFile(f, strng):
    f.write(str(strng))
    return 1
        
def buildNet():
    net = parseInput()
    net.showStartConditions()
    
def parseInput():
    """generate list to populate network class object"""
    f = readFile()
    
    metVec = getMetVec(f)
    netList = getNetList(f)
    runCom = getRunCom(f)
    vEnz = matrix(getConcs(netList))
    vMet = matrix(getConcs(metVec))
    mKConsts = matrix(buildKinConstMat(netList))
    #print(vEnz,vMet,mKConsts)
    nodes = getNodeArray(netList)
    #print(nodes)
    mStoch = getStochMat(nodes)
    #mStoch.display()
    runCom = parseRunCom(getRunCom(f))

    pathway = network([vEnz, vMet, mStoch, mKConsts], runCom)
    

    return pathway
    
def readFile():
    #print("Gimme the filepath")
    #path = userInput
    path = '/Volumes/QUINTIN/Project/MetabolSiM/netList221mM.txt'
    f = open(path, 'r')
    
    red = f.read(536870912) #read in as string, truncate at 1/2GB (should be plenty)
    f.close()
    return red

def getMetVec(red):
    i = 0
    while (i < len(red)) and (red[i:i+4] != '*Met') and (red[i:i+4] != '*met'):
        i += 1
    k = i + 1
    while (k < len(red)) and (red[k] != '*'):
        k += 1
    #print(i,k)
    return red[(i+4):(k-1)]

def getNetList(red):
    i = 0
    while (i < len(red)) and (red[i:i+4] != '*Net') and (red[i:i+4] != '*net'):
        i += 1
    k = i + 1
    while (k < len(red)) and (red[k] != '*'):
        k += 1
    #print(i,k)
    return red[(i+4):(k-1)]

def parseRunCom(strnTemp):
    k = 0
    runCom = []
    while (strnTemp[k] != '(') and (k < len(strnTemp)):          
        k += 1
    i = 0
    while (strnTemp[i] != ')') and (i < len(strnTemp)):
        i += 1
    strng = strnTemp[k+1:i]
    print(strng)
    
    #return strng
    strng = strng.replace(',','\n')
    strng = strng.splitlines()
    print(strng)
    
    if (strng[0] == 't'):
        runCom = [0,float(strng[1]),float(strng[2])]
    else:
        print('boogedy')

    print(runCom)
    return runCom

    

def getRunCom(red):
    i = 0
    while (i < len(red)) and (red[i:i+4] != '*Run') and (red[i:i+4] != '*run'):
        i += 1
    k = i + 1
    while (k < len(red)) and (red[k] != '*'):
        k += 1
    #print(i,k)
    return red[(i+4):(k)]


def getConcs(strng):
    strnTemp = strng.splitlines()
    concs = []
    k = 0
    #for k in strnTemp:
    while (k < len(strnTemp)):
        if '(' in strnTemp[k]:
            concs.append([getConc(strnTemp[k])])
        k+=1

        
    return concs

def getConc(strng):
    i = 0
    while (strng[i] != '('):
        i+=1
        if (i >= len(strng)):
            #There are no user-defined constants defined by bracketing ( [ .. ])
            #Return generic kinetic constants
            conc = 0.05
            return conc            
    k = i
    while (strng[k] != ')') and (k < len(strng)):
        k+=1
    conc = (strng[(i+1):(k)])
    return float(conc)
            
    


def buildKinConstMat(strng):
    strnTemp = strng.splitlines()
    kConsts = []
    k = 0
    #for k in strnTemp:
    while (k < len(strnTemp)):
        if '[' in strnTemp[k]:
            kConsts.append(getKonsts(strnTemp[k]))
        elif '>' in strnTemp[k]:
            print("using default kinetic constants")
            #There are no user-defined constants defined by bracketing ( [ .. ])
            #Return generic kinetic constants
            kConsts.append([.75, 300,1.0, 100])
        k+=1
    return kConsts

def getKonsts(strng):
    i = 0
    while (strng[i] != '['):
        i+=1
        """if (i >= len(strng)):
            print("using default kinetic constants")
            #There are no user-defined constants defined by bracketing ( [ .. ])
            #Return generic kinetic constants
            consts = [.75, 300,1.0, 100]
            return consts  """          
    k = i
    while (strng[k] != ']') and (k < len(strng)):
        k+=1
    #konsts = ((strng[(i+1):(k)]).replace(',','\n')).splitlines()
    konsts = ((strng[(i+1):(k)]).replace(',',' ')).replace(' ','\n').splitlines()
    c = 0
    consts = []
    while (c < len(konsts)):
        if konsts[c]:
            consts.append(float(konsts[c]))
        c+=1
    #while (c < 16):
    #    consts.append(0.01) #fudge factor adds default km's for substrates if not explicitly defined in netlist. maybe replace with error message in final version
    #    c += 1
    return consts

def getNodeArray(f):
    i = 0
    nds = []            #initialize nodelist
    f = f.splitlines()
    while (i<len(f)):
        if (f[i] != '') and ('>' in f[i]):
            nds.append(getNodes(f[i]))  #populate
        i+=1
    #print(nds)                      #debug
    return nds

def getNodes(strng):
    i = 0
    k = 0
    if ('>' in strng):
        while (strng[i] != '>') and (i < len(strng)):
            i += 1
        while ((strng[k] != '[') and (strng[k] != '(') and ( k < (len(strng)-1))):
            k += 1
    
        l = strng[(i+2):k].replace(' ','\n').splitlines()
        #maybe clean up node string here by removing adjacent \n         
        i = 0
        nodes = []
        print(l)
    
        while (i < len(l)):
            nodes.append(int(l[i]))
            i+=1
        return nodes
    
            
def getStochMat(nds):
    i = 0               #initialize i as 0
    m = 0
    n = len(nds)        #number of columns corresponds to number of rxns, e.g. things connected to nodes
    while (i < n):
        m = max(m,max(nds[i]))
        i += 1
    m += 1
    #print (m-1 , "X", n-1)    
    i = 0
    mt = []
    while (i <= m - 1):
        k = 0
        mt.append([])
        while (k < n):
            mt[i].append(0)
            k+=1
        i+=1
    #print(mt)
    i = 0
    #k = 0           #reinit counter
    while (i < n):
        for k in nds[i]:
            #print(k)
            if (k < 0):
                mt[abs(k)][i] += -1
            else:
                mt[k][i] += 1
        #print(i)    
        i += 1
        #print("i = ", i)
    #print(mt,"\n","\n")
    mt.pop(0)       #Remove from list position 0 (this is a row of 0's 
    #global A
    A = matrix(mt)
    #A.display()
    return A

def go():
    A = parseInput(); file = openResFile(); A.tRun(file); file.close()
