# Crystal Plasticity Code.

import math,sys,smallmatrix



########## Loading gauss points for 8 integration points per element ############
def pgauss():
    g = 1.0/math.sqrt(3.0)
    ls = [-1,1,1,-1,-1,1,1,-1]
    lt = [-1,-1,1,1,-1,-1,1,1]
    lz = [-1,-1,-1,-1,1,1,1,1]
    sg = [] ; tg = [] ; zg = [] ; wg = [] 
    for i in range(8):
        s = g*ls[i] ; sg.append(s)
        t = g*lt[i] ; tg.append(t)
        z = g*lz[i] ; zg.append(z)
        w = 1.0     ; wg.append(w)
    return sg,tg,zg,wg
##################### End of loading gauss points ################################

##### Get global shape functions and their derivatives at a given gauss point#####
############################ N_i(X,t) and gradN_i(X,t) ###########################
def SHP3D(ss,tt,zz,XL=[[]]):
    si = [-1,1,1,-1,-1,1,1,-1]
    ti = [-1,-1,1,1,-1,-1,1,1]
    zi = [-1,-1,-1,-1,1,1,1,1]
    SHPL = [[0. for ii in range(8)] for jj in range(4)]
#----------- Local shape functions with derivatives with repect to si,ti and zi
    for i in range(8):
    	N1 = 1.0+ss*si[i]
        N2 = 1.0+tt*ti[i]
        N3 = 1.0+zz*zi[i]
        SHPL[3][i] = N1*N2*N3/8.0
        SHPL[0][i] = si[i]*N2*N3/8.0
        SHPL[1][i] = ti[i]*N1*N3/8.0
        SHPL[2][i] = zi[i]*N1*N2/8.0
#--------- Get global values in refrence configuration

    xsi = smallmatrix.SmallMatrix(3,3)
    isi = smallmatrix.SmallMatrix(3,3)

    xsi.clear()
    isi.clear()

    for i in range(3):
        for j in range(3):
            isi[i,j] = 0.0
    isi[0,0]=isi[1,1]=isi[2,2]=1.0


    XS = [[0. for ii in range(3)] for jj in range(3)]
    SX = [[0. for ii in range(3)] for jj in range(3)]
 

    for i in range(3):
        for j in range(3):
            for k in range(8):
            	xsi[i,j]+=XL[k][j]*SHPL[i][k]

    xsi.solve(isi)
    for i in range(3):
       for j in range(3):
           SX[i][j] = isi[i,j]
#           XS[i][j] = xsi[i,j]
    
    SHPG = [[0. for ii in range(8)] for jj in range(3)]
    for i in range(8):
        T0 = SX[0][0]*SHPL[0][i]+SX[0][1]*SHPL[1][i]+SX[0][2]*SHPL[2][i]
        T1 = SX[1][0]*SHPL[0][i]+SX[1][1]*SHPL[1][i]+SX[1][2]*SHPL[2][i]
        T2 = SX[2][0]*SHPL[0][i]+SX[2][1]*SHPL[1][i]+SX[2][2]*SHPL[2][i]
        SHPG[0][i] = T0
        SHPG[1][i] = T1
        SHPG[2][i] = T2
                
    return SHPG
############### End of calculation of shape functions and drivatives ##################

########### Calculate deformation gradient (dx/dX=sum(x*gradN_i(X,t))) ################
def calc_F(cord,SHP):
    
    dxt_dx0 = 0.0 ; dxt_dy0 = 0.0 ; dxt_dz0 = 0.0
    dyt_dx0 = 0.0 ; dyt_dy0 = 0.0 ; dyt_dz0 = 0.0
    dzt_dx0 = 0.0 ; dzt_dy0 = 0.0 ; dzt_dz0 = 0.0

    for I in range(8):
	dxt_dx0 += SHP[0][I]*cord[I][0]
	dxt_dy0 += SHP[1][I]*cord[I][0]
	dxt_dz0 += SHP[2][I]*cord[I][0]
	
	dyt_dx0 += SHP[0][I]*cord[I][1]
	dyt_dy0 += SHP[1][I]*cord[I][1]
	dyt_dz0 += SHP[2][I]*cord[I][1]

	dzt_dx0 += SHP[0][I]*cord[I][2]
	dzt_dy0 += SHP[1][I]*cord[I][2]
	dzt_dz0 += SHP[2][I]*cord[I][2]
 
    FG = [[0. for ii in range(3)] for jj in range(3)]
    FG[0][0] = dxt_dx0 ; FG[0][1] = dxt_dy0 ; FG[0][2] = dxt_dz0
    FG[1][0] = dyt_dx0 ; FG[1][1] = dyt_dy0 ; FG[1][2] = dyt_dz0
    FG[2][0] = dzt_dx0 ; FG[2][1] = dzt_dy0 ; FG[2][2] = dzt_dz0

    return FG
#################### End of calculation of deformation gradient ###################

############ Calculate A(n,n)*B(n,n)############
def matmul(len,A=[[]],B=[[]]):
#    BT = zip(*B)
    C = [[0. for ii in range(len)] for jj in range(len)]
    for i in range(len):
        for j in range(len):
            for k in range(len):
                C[i][j] += A[i][k]*B[k][j]
    return C

###########################################################################
##################### Schmid Tensor for FCC materials #####################
###########################################################################

def calc_schmid(crystal_type,n_slip,qrot):
    schmid = [[[0. for i in range(3)] for j in range(3)] for k in range(n_slip)]

    r1=1.0
    r2=1.0/math.sqrt(2.0)
    r3=1.0/math.sqrt(3.0)
    r6=1.0/math.sqrt(6.0)
 
    if crystal_type is 'fcc':

        cns = [[0. for ii in range(3)] for jj in range(n_slip)]
        cms = [[0. for ii in range(3)] for jj in range(n_slip)]
        vec1 = [[0. for ii in range(3)] for jj in range(1)]
        vec2 = [[0. for ii in range(1)] for jj in range(3)]

#---------------fcc slip system 1----------
	ind_slip = 0

	cns[ind_slip][0] = r3
	cns[ind_slip][1] = r3
	cns[ind_slip][2] = r3

        cms[ind_slip][0] = 0.0
        cms[ind_slip][1] = r2
        cms[ind_slip][2] = -r2

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
	

	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]

#---------------fcc slip system 2----------
	ind_slip = 1

	cns[ind_slip][0] = r3
	cns[ind_slip][1] = r3
	cns[ind_slip][2] = r3

        cms[ind_slip][0] = -r2
        cms[ind_slip][1] = 0.0
        cms[ind_slip][2] = r2

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
	

	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]


        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]
                     
#---------------fcc slip system 3----------
	ind_slip = 2

	cns[ind_slip][0] = r3
	cns[ind_slip][1] = r3
	cns[ind_slip][2] = r3

        cms[ind_slip][0] = r2
        cms[ind_slip][1] = -r2
        cms[ind_slip][2] = 0.0

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
	

	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]
                     
#---------------fcc slip system 4----------
	ind_slip = 3

	cns[ind_slip][0] = r3
	cns[ind_slip][1] = -r3
	cns[ind_slip][2] = -r3

        cms[ind_slip][0] = 0.0
        cms[ind_slip][1] = -r2
        cms[ind_slip][2] = r2

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
	

	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]
                     
#---------------fcc slip system 5----------
	ind_slip = 4

	cns[ind_slip][0] = r3
	cns[ind_slip][1] = -r3
	cns[ind_slip][2] = -r3

        cms[ind_slip][0] = -r2
        cms[ind_slip][1] = 0.0
        cms[ind_slip][2] = -r2

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
	

	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]

#---------------fcc slip system 6----------
	ind_slip = 5

	cns[ind_slip][0] = r3
	cns[ind_slip][1] = -r3
	cns[ind_slip][2] = -r3

        cms[ind_slip][0] = r2
        cms[ind_slip][1] = r2
        cms[ind_slip][2] = 0.0

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
	

	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
        
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]

#---------------fcc slip system 7----------
	ind_slip = 6

	cns[ind_slip][0] = -r3
	cns[ind_slip][1] = r3
	cns[ind_slip][2] = -r3

        cms[ind_slip][0] = 0.0
        cms[ind_slip][1] = r2
        cms[ind_slip][2] = r2

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
    
           
	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]
                     
#---------------fcc slip system 8----------
	ind_slip = 7

	cns[ind_slip][0] = -r3
	cns[ind_slip][1] = r3
	cns[ind_slip][2] = -r3

        cms[ind_slip][0] = r2
        cms[ind_slip][1] = 0.0
        cms[ind_slip][2] = -r2

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
	

	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]
 
#---------------fcc slip system 9----------
	ind_slip = 8

	cns[ind_slip][0] = -r3
	cns[ind_slip][1] = r3
	cns[ind_slip][2] = -r3

        cms[ind_slip][0] = -r2
        cms[ind_slip][1] = -r2
        cms[ind_slip][2] = 0.0

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
	

	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]

#---------------fcc slip system 10----------
	ind_slip = 9

	cns[ind_slip][0] = -r3
	cns[ind_slip][1] = -r3
	cns[ind_slip][2] = r3

        cms[ind_slip][0] = 0.0
        cms[ind_slip][1] = -r2
        cms[ind_slip][2] = -r2

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
	

	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]

#---------------fcc slip system 11----------
	ind_slip = 10

	cns[ind_slip][0] = -r3
	cns[ind_slip][1] = -r3
	cns[ind_slip][2] = r3

        cms[ind_slip][0] = r2
        cms[ind_slip][1] = 0.0
        cms[ind_slip][2] = r2

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
	

	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]

#---------------fcc slip system 12----------
	ind_slip = 11

	cns[ind_slip][0] = -r3
	cns[ind_slip][1] = -r3
	cns[ind_slip][2] = r3

        cms[ind_slip][0] = -r2
        cms[ind_slip][1] = r2
        cms[ind_slip][2] = 0.0

	for ij in range(3):
           vec2[ij][0] = cns[ind_slip][ij]
	

	for ik in range(3):
            vec1[0][ik] = cms[ind_slip][ik]
	
        tmp_v1 = [[0. for ii in range(3)] for jj in range(1)]
	tmp_v2 = [[0. for ii in range(1)] for jj in range(3)]
        tmp1 = [[0. for ii in range(3)] for jj in range(3)]

        for i in range(3):
            for j in range(3):
                tmp_v1[0][i] += qrot[i][j]*vec1[0][j]
                tmp_v2[i][0] += qrot[i][j]*vec2[j][0]

        for i in range(3):
            for j in range(3):
                tmp1[i][j] = tmp_v1[0][i]*tmp_v2[j][0]
            
        for i in range(3):
            for j in range(3):
                schmid[ind_slip][i][j] = tmp1[i][j]
    
    return schmid

############## End of Schmid Tensor for FCC materials #####################


############## Rotation Matrix based on Euler angles #####################

def euler_rot(angle,qrot):
    const_pi = math.acos(-1.0)
  
    phi = angle[0]*const_pi/180.0
    theta = angle[1]*const_pi/180.0
    omega = angle[2]*const_pi/180.0
	
    sp = math.sin(phi)                      
    cp = math.cos(phi)                     
    st = math.sin(theta)                     
    ct = math.cos(theta)                    
    so = math.sin(omega)                    
    co = math.cos(omega)   

    qrot[0][0] = co*cp-so*sp*ct
    qrot[1][0] = co*sp+so*ct*cp   
    qrot[2][0] = so*st   
    qrot[0][1] = -so*cp-sp*co*ct 
    qrot[1][1] = -so*sp+ct*co*cp
    qrot[2][1] = co*st
    qrot[0][2] = sp*st       
    qrot[1][2] = -st*cp       
    qrot[2][2] = ct


    return qrot

###################### End of Rotation Matrix #####################

############## Elastic Material 2nd and 4th order tensor ###################

def cmat(c11,c12,c44,qrot):

    cij = [[0. for i in range(6)] for j in range(6)] 
    
#---------------- cij for cubic materials 
    cij[0][0] = cij[1][1] = cij[2][2] = c11
    cij[0][1] = cij[0][2] = cij[1][0] = cij[1][2] = cij[2][0] = cij[2][1] = c12
    cij[3][3] = cij[4][4] = cij[5][5] = c44


    cijkl = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] 

#----------------- Get cijkl  
    for i in range(3):
        for j in range(3):
            cijkl[i][i][j][j] = cij[i][j]

    for i in range(3):
        cijkl[i][i][0][1] = cij[i][3]
        cijkl[i][i][1][0] = cij[i][3]
        cijkl[i][i][1][2] = cij[i][4]
        cijkl[i][i][2][1] = cij[i][4]
        cijkl[i][i][0][2] = cij[i][5]
        cijkl[i][i][2][0] = cij[i][5]
        
        cijkl[0][1][i][i] = cijkl[1][0][i][i] = cij[3][i]
        cijkl[1][2][i][i] = cijkl[2][1][i][i] = cij[4][i]
        cijkl[0][2][i][i] = cijkl[2][0][i][i] = cij[5][i]
    
    cijkl[0][1][0][1] = cijkl[0][1][1][0] = cijkl[1][0][0][1] = cijkl[1][0][1][0] = cij[3][3]
    cijkl[0][1][1][2] = cijkl[0][1][2][1] = cijkl[1][0][1][2] = cijkl[1][0][2][1] = cij[3][4]
    cijkl[0][1][0][2] = cijkl[0][1][2][0] = cijkl[1][0][0][2] = cijkl[1][0][2][0] = cij[3][5]
    cijkl[1][2][0][1] = cijkl[1][2][1][0] = cijkl[2][1][0][1] = cijkl[2][1][1][0] = cij[4][3]
    cijkl[1][2][1][2] = cijkl[1][2][2][1] = cijkl[2][1][1][2] = cijkl[2][1][2][1] = cij[4][4]
    cijkl[1][2][0][2] = cijkl[1][2][2][0] = cijkl[2][1][0][2] = cijkl[2][1][2][0] = cij[4][5]
    cijkl[0][2][0][1] = cijkl[0][2][1][0] = cijkl[2][0][0][1] = cijkl[2][0][1][0] = cij[5][3]
    cijkl[0][2][1][2] = cijkl[0][2][2][1] = cijkl[2][0][1][2] = cijkl[2][0][2][1] = cij[5][4]
    cijkl[0][2][0][2] = cijkl[0][2][2][0] = cijkl[2][0][0][2] = cijkl[2][0][2][0] = cij[5][5]


#----------------- C_MAT for Crytal with respect to orientation 
    C_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] 
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    
                    for m in range(3):
                        for n in range(3):
                            for p in range(3):
                                for q in range(3):
                                    C_mat[i][j][k][l] += qrot[i][m]*qrot[j][n]*qrot[k][p]*qrot[l][q]*cijkl[m][n][p][q]


    return C_mat
########################## End of calculation C_mat ##########################

########################## Calculate matrix [A] ##############################
def calc_A(F,Fp):
    F_tau = [[0. for ii in range(3)] for jj in range(3)]
    FT_tau = [[0. for ii in range(3)] for jj in range(3)]
    FpI_t = [[0. for ii in range(3)] for jj in range(3)]
    FpIT_t = [[0. for ii in range(3)] for jj in range(3)]
    TMP1 = [[0. for ii in range(3)] for jj in range(3)]
    TMP2 = [[0. for ii in range(3)] for jj in range(3)]
    A_mat = [[0. for ii in range(3)] for jj in range(3)]
    Fp_t = smallmatrix.SmallMatrix(3,3)
    Ip = smallmatrix.SmallMatrix(3,3)
    Fp_t.clear()
    Ip.clear()
    for i in range(3):
        for j in range(3):
            Ip[i,j] = 0.0
            Ip[0,0]=Ip[1,1]=Ip[2,2]=1.0


    for i in range(3):
        for j in range(3):
            Fp_t[i,j] = Fp[i][j]

    Fp_t.solve(Ip)

    for i in range(3):
        for j in range(3):
            FpI_t[i][j] = Ip[i,j]
            FpIT_t[i][j] = Ip[j,i]


    for i in range(3):
        for j in range(3):
            F_tau[i][j] = F[i][j]
            FT_tau[i][j] = F[j][i]

    TMP1 = matmul(3,F_tau,FpI_t)
    TMP2 = matmul(3,FT_tau,TMP1)
    A_mat = matmul(3,FpIT_t,TMP2)

    return A_mat
########################## End of A_mat ##################################

################ Calculate matrix [B_alpha] ##############################
def calc_B(n_slip,A_mat,schmid):
    B_alpha = [[[0. for i in range(3)] for j in range(3)] for k in range(n_slip)]
    TMP1 = [[0. for ii in range(3)] for jj in range(3)]
    TMP1T = [[0. for ii in range(3)] for jj in range(3)]
    TMP20 = [[0. for ii in range(3)] for jj in range(3)]
    TMP21 = [[0. for ii in range(3)] for jj in range(3)]
    TMP2 = [[0. for ii in range(3)] for jj in range(3)]
    for k in range(n_slip):
        for j in range(3):
            for i in range(3):
                TMP1[i][j] = schmid[k][i][j]
                TMP1T[i][j] = schmid[k][j][i]
                

        TMP20 = matmul(3,A_mat,TMP1)
        TMP21= matmul(3,TMP1T,A_mat)

        for it in range(3):
            for jt in range(3):
                TMP2[it][jt] = TMP20[it][jt]+TMP21[it][jt]
        
        for i in range(3):
            for j in range(3):
                B_alpha[k][i][j]=TMP2[i][j]
        
    return B_alpha
###################### End of B_alpha #####################################
    
############# Calculate the trial elastic stress S_trial ##################
def stress_trial(C_mat,A_mat):
    delta_kron = [[0. for ii in range(3)] for jj in range(3)]
    S_trial = [[0. for ii in range(3)] for jj in range(3)]
    
    delta_kron[0][0] = delta_kron[1][1]= delta_kron[2][2] = 1.0
	
    for i in range(3):
        for j in range(3):
            for m in range(3):
                for n in range(3):
                    S_trial[i][j] += C_mat[i][j][m][n]*0.5*(A_mat[m][n]-delta_kron[m][n])

    return S_trial
########################## End of S_trial ####################################

################ Calculate the matrix C_ALPHA for each slip system ###########
def calc_C_alpha(n_slip,C_mat,B_alpha):
    C_alpha = [[[0. for i in range(3)] for j in range(3)] for k in range(n_slip)]    

    for k in range(n_slip):
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    for n in range(3):
                        C_alpha[k][i][j] += C_mat[i][j][m][n]*0.5*B_alpha[k][m][n]

    return C_alpha


####################### End of C_alpha ####################################

################  Resolved Shear Stress ##########
def resolvedshear(n_slip,stress,schmidfactor):
    tau_alpha = [0.0 for k in range(n_slip)]
    for k in range(n_slip):
        for i in range(3):
            for j in range(3):
                tau_alpha[k] += stress[i][j]*schmidfactor[k][i][j]
 
    return tau_alpha

#################### End of Resolved Shear Stress ##############

#################  Reduce 4th order 3*3*3*3 tensor to 2nd 6*6 matrix ########

def reduce_mat(as_mat): 
       as_redu = [[0.0 for i in range(6)] for j in range(6)]

       for i in range(3):
           for j in range(3):
               as_redu[i][j] = as_mat[i][i][j][j]

       for i in range(3):
           as_redu[i][3] = as_mat[i][i][0][1]+as_mat[i][i][1][0]
           as_redu[i][4] = as_mat[i][i][1][2]+as_mat[i][i][2][1]
           as_redu[i][5] = as_mat[i][i][0][2]+as_mat[i][i][2][0]

       for j in range(3):
           as_redu[3][j] = as_mat[0][1][j][j]
           as_redu[4][j] = as_mat[1][2][j][j]
           as_redu[5][j] = as_mat[0][2][j][j]

       as_redu[3][3] = as_mat[0][1][0][1]+as_mat[0][1][1][0]
       as_redu[3][4] = as_mat[0][1][1][2]+as_mat[0][1][2][1]
       as_redu[3][5] = as_mat[0][1][0][2]+as_mat[0][1][2][0]

       as_redu[4][3] = as_mat[1][2][0][1]+as_mat[1][2][1][0]
       as_redu[4][4] = as_mat[1][2][1][2]+as_mat[1][2][2][1]
       as_redu[4][5] = as_mat[1][2][0][2]+as_mat[1][2][2][0]

       as_redu[5][3] = as_mat[0][2][0][1]+as_mat[0][2][1][0]
       as_redu[5][4] = as_mat[0][2][1][2]+as_mat[0][2][2][1]
       as_redu[5][5] = as_mat[0][2][0][2]+as_mat[0][2][2][0]

       return as_redu

    
#################  end of Reduce 4th order 3*3*3*3 tensor to 2nd 6*6 matrix ########

#################  Reduce 2nd order 3*3 tensor to 1st 6 vector ########

def reduce_vec(as_star):
    as_vec = [0.0 for i in range(6)]
    for i in range(3):
        as_vec[i] = as_star[i][i]
    
    as_vec[3] = 0.5*(as_star[0][1]+as_star[1][0])
    as_vec[4] = 0.5*(as_star[1][2]+as_star[2][1])
    as_vec[5] = 0.5*(as_star[0][2]+as_star[2][0])

    return as_vec

#################  end of Reduce 2nd order 3*3 tensor to 1st 6 vector ########


######### Inflate 1st 6 vector to 2nd 3*3 tensor ################
def inflate_vec(as_vec):
    as_star = [[0.0 for i in range(3)] for j in range(3)]
    for i in range(3):
        as_star[i][i] = as_vec[i]

    as_star[0][1] = as_vec[3]
    as_star[1][2] = as_vec[4]    
    as_star[0][2] = as_vec[5]

    as_star[1][0] = as_vec[3]
    as_star[2][1] = as_vec[4]
    as_star[2][0] = as_vec[5]

    return as_star

######### end of Inflate 1st 6 vector to 2nd 3*3 tensor ################


################# Calculation of GN vector ###################

def calc_GNvec(n_slip,S_star,S_trial,dgam,C_alpha):
    GN_vec = [0.0 for i in range(6)]
    TMP1 = [[0.0 for i in range(3)] for j in range(3)]

    for i  in range(3):
        for j in range(3):
                TMP1[i][j] = 0.0
                for k in range(n_slip):
                    TMP1[i][j] += dgam[k]*C_alpha[k][i][j]

    for i in range(3):
        for j in range(3):
            TMP1[i][j] += S_star[i][j]-S_trial[i][j]

    GN_vec = reduce_vec(TMP1)

    return GN_vec
    

################# end of Calculation of GN vector ###################



class Node:
    def __init__(self,inode,coorx,coory,coorz):
        self.inode = inode
        self.coorx = coorx
        self.coory = coory
        self.coorz = coorz
    def __repr__(self):
        return "(%d,%f,%f,%f)" % (self.inode,self.coorx,self.coory,self.coorz)

class Element:
    def __init__(self,idx,nodelist=[]):
        self.index = idx
        self.nodes = nodelist[:]

    def __repr__(self):
        return "Element(%d,%s)" % (self.index, self.nodes)


####################### Flow rule #########################
class update_statev:
    def __init__(self,tau_alpha,res_ssd):
        self.tau_alpha = tau_alpha
        self.res_ssd = res_ssd

    def RES(self):


        const_w1 = 1.0 ; const_w2 = 1.4
        const_ss = 377100000.0
        const_a = 1.115
        const_h0 = 1501000000.0
        const_m = 0.01
        dt = 0.1
        g0dot = 24530000.0

        res0 = [0.0 for i in range(len(self.tau_alpha))]
        res_t = [0.0 for i in range(len(self.tau_alpha))]
        dgamma = [0.0 for i in range(len(self.tau_alpha))]
        dgam_dtau = [0.0 for i in range(len(self.tau_alpha))]

        for i in range(len(self.tau_alpha)):
            res0[i] = self.res_ssd[i]
#            dgamma[i] = 1.0e-4
    

        for k in range(len(self.tau_alpha)):
            for i in range(len(self.tau_alpha)):
                if i == k:
                    const_qab = const_w1
                elif i != k:
                    const_qab = const_w2
                
                ratio_res=1.0-(res0[i]/const_ss)

                res_t[k] += const_qab*const_h0*abs(dgamma[i])*(ratio_res)**const_a


        for k in range(len(self.tau_alpha)):
            self.res_ssd[k] = res_t[k]+res0[k]


        for k in range(len(self.tau_alpha)):
            if self.res_ssd[k] >= 1.0:
                ratio_alpha = self.tau_alpha[k]/self.res_ssd[k]
                
            
                if self.tau_alpha[k] >= 0.0:
                    const_sign = 1.0
                elif self.tau_alpha[k] < 0.0:
                    const_sign = -1.0
	
                m_inv = 1.0/const_m
                dgamma[k] = dt*const_sign*g0dot*(abs(ratio_alpha))**m_inv
                
                

######	Calculation of dgamma_dtau (The sgn(Tau_alpha) is multiplied with another sgn(Tau_alpha) from the derivative)
                res_inv = 1.0/self.res_ssd[k]
                dgam_dtau[k] = dt*res_inv*g0dot*m_inv*(abs(ratio_alpha))**(m_inv-1.0)

            elif self.res_ssd[k] < 1.0:
                dgamma[k] = 0.0
                dgam_dtau[k] = 0.0

    
	
        return dgamma,dgam_dtau,self.res_ssd

################ User Material ###########

class UMAT:
    def __init__(self,S_trial,C_mat,Fp,C_alpha,schmid,F_t,F_tau,res):
        self.S_trial = S_trial
        self.C_maT = C_mat
        self.Fp = Fp
        self.C_alpha = C_alpha
        self.schmid = schmid
        self.F_t = F_t
        self.F_tau = F_tau
        self.res = res
#        self.tau_alpha = [0.0 for k in range(len(self.schmid))]

#---------------------- Polar Decomposition of Rlative Gradient Deformation -------------------
    def polardecomp(self):

#-------------- Calculation of relative deformation gradient F_REL = F_tau*(F_t)^-1 -----------
        F_REL = [[0.0 for i in range(3)] for j in range(3)]
        TMPrel = smallmatrix.SmallMatrix(3,3)
        Irel = smallmatrix.SmallMatrix(3,3)
        TMPrel.clear()
        Irel.clear()
        for i in range(3):
            for j in range(3):
                Irel[i,j] = 0.0
        Irel[0,0]=Irel[1,1]=Irel[2,2]=1.0

        for i in range(3):
            for j in range(3):
                TMPrel[i,j] = self.F_t[i][j]

        TMPrel.solve(Irel)

        for i in range(3):
            for j in range(3):
                F_REL[i][j] = 0.0
                for k in range(3):
                    F_REL[i][j] += self.F_tau[i][k]*Irel[k,j]
                
#----------------------- square root of a positive matrix U=sqrt(C) ------------

        R = [[0. for i in range(3)] for j in range(3)]
        C = [[0. for i in range(3)] for j in range(3)]
        Csquare = [[0. for i in range(3)] for j in range(3)]
        Iden = [[0. for i in range(3)] for j in range(3)]
        U = [[0. for i in range(3)] for j in range(3)]
        invU = [[0. for i in range(3)] for j in range(3)]

        Iden[0][0] = Iden[1][1] = Iden[2][2] = 1.0

#----------------------- C = FTF --------------------------------------
        for i in range(3):
                for j in range(3):
                    for k in range(3):
                        C[i][j] += F_REL[k][i]*F_REL[k][j]

#---------------------- C^2 ------------------------------------------
        for i in range(3):
                for j in range(3):
                    for k in range(3):
                        Csquare[i][j] += C[i][k]*C[k][j]


        I_C = C[0][0]+C[1][1]+C[2][2]

#--------------------- Invarients of C -------------------------------
        I_Csquare = Csquare[0][0]+Csquare[1][1]+Csquare[2][2] 
        II_C = 0.5*(I_C**2-I_Csquare)
        III_C = C[0][0]*(C[1][1]*C[2][2]-C[1][2]*C[2][1])-\
            C[0][1]*(C[1][0]*C[2][2]-C[2][0]*C[1][2])+\
            C[0][2]*(C[1][0]*C[2][1]-C[1][1]*C[2][0])

        k = I_C**2-3.0*II_C

        if k<0.0001:
            lamda = math.sqrt(I_C/3)
            for i in range(3):
                U[i][i] = lamda
                invU[i][i] = 1/lamda

        elif k>=0.0001:

            l = I_C*(I_C**2-4.5*II_C)+13.5*III_C

            teta = (1.0/3.0)*math.acos(l/(k**(1.5)))

            alpha = 2.0*math.pi/3.0

    #--------------------- Eigenvalues of U ----------------------------- 
            lamda1 = math.sqrt((1.0/3.0)*(I_C+2.0*math.sqrt(k)*math.cos(teta)))
            lamda2 = math.sqrt((1.0/3.0)*(I_C+2.0*math.sqrt(k)*math.cos(alpha+teta)))
            lamda3 = math.sqrt((1.0/3.0)*(I_C+2.0*math.sqrt(k)*math.cos(2.0*alpha+teta)))


            #III_U = math.sqrt(III_C)
            #I_U = lamda+math.sqrt(-(lamda**2)+I_C+2.0*III_C/lamda)
            #II_U = (I_U**2-I_C)/2.0

    #--------------------- Invarients of U ---------------------------
            I_U = lamda1+lamda2+lamda3
            II_U = lamda1*lamda2+lamda1*lamda3+lamda2*lamda3
            III_U = lamda1*lamda2*lamda3

    #--------------------- U and Inverse of U ------------------------

            for i in range(3):
                    for j in range(3):
                            U[i][j] = (1.0/(I_U*II_U-III_U))*(I_U*III_U*Iden[i][j]+\
                                                                      (I_U**2-II_U)*C[i][j]-Csquare[i][j])
                            invU[i][j] = (1.0/III_U)*(II_U*Iden[i][j]-\
                                                                      I_U*U[i][j]+C[i][j])
#----------------------- R = FU^-1 ----------------------------- 
        for i in range(3):
                for j in range(3):
                    for k in range(3):
                        R[i][j] += F_REL[i][k]*invU[k][j]

#---------------------- End of Polar Decomposition ---------------

        return U,R

#----------------------- Storing Fp into Fp_t-------------------
    def fpt(self):
        Fp_t = [[0.0 for i in range(3)] for j in range(3)]
        Fp_t = self.Fp
        return Fp_t
#--------------------------------------------------------------

#------------ Forth order Kronecker Delta ---------------------
    def delta_kron4(self):
        delta_kron = [[0. for ii in range(3)] for jj in range(3)]
        delta_kron[0][0] = delta_kron[1][1]= delta_kron[2][2] = 1.0

        delta_kron4d = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)] 

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        delta_kron4d[i][j][k][l] += delta_kron[i][k]*delta_kron[j][l]+delta_kron[i][l]*delta_kron[j][k]

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        delta_kron4d[i][j][k][l] = 0.5*delta_kron4d[i][j][k][l]
        
        return delta_kron4d
#---------------------------------------------------------------

#------------------ Iteration process for time integration of Fp --------        

    def itr(self):
        n_slip = len(self.schmid)

        tau_alpha = [0.0 for k in range(n_slip)]
        tau_alpha = resolvedshear(n_slip,self.S_trial,self.schmid)

        res_ssd = self.res
        
        
#---------------------- Update state variables ----------------
#????????????????? Get and store previous values for the slip system resistances????????
        
        s1 = update_statev(tau_alpha,res_ssd)
        u1 = s1.RES()
        dgam = u1[0]
        dgma_dta = u1[1]
        res_ssd = u1[2]
        

        TMP1 = [[0.0 for it in range(3)] for jt in range(3)]
        for it in range(3):
            for jt in range(3):
                for kt in range(n_slip):
                    TMP1[it][jt] += dgam[kt]*self.C_alpha[kt][it][jt]


#-------------------- Calculation of first value of 2nd PKS-----------------
        S_star = [[0.0 for it in range(3)] for jt in range(3)]
        for it in range(3):
            for jt in range(3):
                S_star[it][jt] = self.S_trial[it][jt]-TMP1[it][jt]

#----------------- iteration loop for S_star -------------------
        niter = 20
        iitr = 0
        ratio_norm = 1.0e10
        while (iitr <= 20 and abs(ratio_norm) >= 0.0001):

 #------------ Get the resolved shear stress and update state variables with S_star-----------------
            tau_alpha = [0.0 for k in range(n_slip)]
            tau_alpha = resolvedshear(n_slip,S_star,self.schmid)


            s2 = update_statev(tau_alpha,res_ssd)
            u2 = s2.RES()
            dgam = u2[0]
            dgam_dta = u2[1]
            res_ssd = u2[2]
            
#---------
            GT_mat = [[[0. for i in range(3)] for j in range(3)] for k in range(n_slip)]

            for k in range(n_slip):
                for i in range(3):
                    for j in range(3):
                        GT_mat[k][i][j] = 0.5*(self.schmid[k][i][j]+self.schmid[k][j][i])*dgam_dta[k]
#-----------
#----------
            RJ_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
            for i in range(3):
                for j in range(3):
                    for m in range(3):
                        for n in range(3):                  
                            for k in range(n_slip):
                                RJ_mat[i][j][m][n] += self.C_alpha[k][i][j]*GT_mat[k][m][n]
            
#-------------
            delta4 = UMAT.delta_kron4(self)
            for i in range(3):
                for j in range(3):
                    for m in range(3):
                        for n in range(3):
                            RJ_mat[i][j][m][n] += delta4[i][j][m][n]

            RJ_reduced = reduce_mat(RJ_mat)
#--------------
#------------

            RJ = smallmatrix.SmallMatrix(6,6)
            IRJ = smallmatrix.SmallMatrix(6,6)
            RJ.clear()
            IRJ.clear()
            for i in range(6):
                for j in range(6):
                    IRJ[i,j] = 0.0
            IRJ[0,0] = IRJ[1,1] = IRJ[2,2] = IRJ[3,3] = IRJ[4,4] = IRJ[5,5] = 1.0


            for i in range(6):
                for j in range(6):
                    RJ[i,j] = RJ_reduced[i][j]

            RJ.solve(IRJ)

            for i in range(6):
                for j in range(6):
                    RJ_reduced[i][j] = IRJ[i,j]
#---------------
#-------------  
            GN_vec = calc_GNvec(n_slip,S_star,self.S_trial,dgam,self.C_alpha)
#----------------
#-----------------

            S_vec = reduce_vec(S_star)

            TMP_RG = [0.0 for i in range(6)]
            TMP_vec = [0.0 for i in range(6)]

            for i in range(6):
                for j in range(6):
                    TMP_RG[i] += RJ_reduced[i][j]*GN_vec[j]

            for i in range(6):
                TMP_vec[i] = S_vec[i]

            for i in range(6):
                S_vec[i] -= TMP_RG[i]
            

            S_star = inflate_vec(S_vec)
            
        
            dot_product_S_vec = 0.0
            for i in range(6):
                dot_product_S_vec += S_vec[i]*S_vec[i]

            rnorm_s_vec = math.sqrt(dot_product_S_vec)


            dot_product_TMP_vec = 0.0
            for i in range(6):
                dot_product_TMP_vec += TMP_vec[i]*TMP_vec[i]

            rnorm_tmp_vec = math.sqrt(dot_product_TMP_vec)


            if abs(rnorm_tmp_vec) > 1.0:
                ratio_norm = (rnorm_s_vec-rnorm_tmp_vec)/rnorm_tmp_vec
#--------------------

            print iitr+1,ratio_norm

            iitr += 1



        return u1
    
#------------------ End of Iteration method for time integration of Fp --------          
#--------------------------------------------------------------





class Mesh:
    def __init__(self,xelements=2,yelements=2,zelements=2):
        self.nodelist = []
        self.ellist = [] 

        dx = 1.0/xelements
        dy = 1.0/yelements
        dz = 1.0/zelements


        node_index = 0
        for i in range(zelements+1):
            for j in range(yelements+1):
                for k in range(xelements+1):
                    node = Node(node_index,k*dx,j*dy,i*dz)
                    self.nodelist.append(node)
                    node_index += 1
                    
        element_index = 0
        for i in range(zelements):
            for j in range(yelements):
                for k in range(xelements):
                    nodes = [self.nodelist[i*(xelements+1)*(yelements+1)+j*(xelements+1)+k],
                             self.nodelist[i*(xelements+1)*(yelements+1)+j*(xelements+1)+k+1],
                             self.nodelist[i*(xelements+1)*(yelements+1)+(j+1)*(xelements+1)+k+1],
                             self.nodelist[i*(xelements+1)*(yelements+1)+(j+1)*(xelements+1)+k],
                             self.nodelist[(i+1)*(xelements+1)*(yelements+1)+j*(xelements+1)+k],
                             self.nodelist[(i+1)*(xelements+1)*(yelements+1)+j*(xelements+1)+k+1],
                             self.nodelist[(i+1)*(xelements+1)*(yelements+1)+(j+1)*(xelements+1)+k+1],
                             self.nodelist[(i+1)*(xelements+1)*(yelements+1)+(j+1)*(xelements+1)+k]]
                            
                    self.ellist.append(Element(element_index,nodes))
                    element_index += 1



        
############ starting element loop #####################
    def solve(self):
        XYZ = [[0.0 for i in range(3)] for j in range(8)]
        xyz = [[0.0 for i in range(3)] for j in range(8)]
        gdisp = [[0.0 for i in range(3)] for j in range(len(self.nodelist))]
        gdisp[26][0] = 0.1
    
        disp = [0.0 for i in range(24)]

        for ielem in self.ellist:
                
            qrot = [[0. for ii in range(3)] for jj in range(3)]

#            qrot_t = [[0. for ii in range(3)] for jj in range(3)]

            qrot[0][0] = 1.0 ; qrot[1][1] = 1.0 ; qrot[2][2] = 1.0

            angle = [0. for ii in range(3)]
            angle[0] = 45.0 ; angle[1] = 30.0 ; angle[2] = 0.0

            qrot = euler_rot(angle,qrot)
 
#            qrot = zip(*qrot_t)
            

#------------ Coordinates in the refrence configuration
            for i in range(1):  ######### 8 to 1
                XYZ[i][0] = ielem.nodes[i].coorx
                XYZ[i][1] = ielem.nodes[i].coory
                XYZ[i][2] = ielem.nodes[i].coorz
                
                for j in range(3):
                    disp[i*3+j] = gdisp[ielem.nodes[i].inode][j]  
                    xyz[i][j] = XYZ[i][j]+gdisp[ielem.nodes[i].inode][j]

            crystal_type = 'fcc' 
            n_slip = 12
            schmid = calc_schmid(crystal_type,n_slip,qrot)


            c11 = 201700000000.000 ; c12 = 134400000000.000 ; c44 = 104500000000.000 
            C_mat = [[[[0. for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]

            C_mat = cmat(c11,c12,c44,qrot)

            gaussw = pgauss()

### ?????????? transfer plastic deformation gradient from state variables to local variables ### 
            
#---------------------- Start Calculation for Gauss Points

#------- Get the Deformation Gradient at time = tau for gauss points for the given element
            
            for ig in range(1): ########## 8 to 1

#???????????????? Transfer deformation gradient to the previous value ?????????

#---------------------- F_t = F_tau previous -------------

                F_t = [[0. for ii in range(3)] for jj in range(3)]
                F_t[0][0] = F_t[1][1] = F_t[2][2] = 1.0

                F_t[0][0] = 0.999966679104601 ;F_t[0][1] = 6.609837380555095e-4 ;F_t[0][2] = -4.687626276934642e-4 
                F_t[1][0] = 4.518037680970892e-7 ;F_t[1][1] = 1.00154297191995 ;F_t[1][2] = 2.184383201542781e-7  
                F_t[2][0] = -3.480313499383450e-6 ;F_t[2][1] = 5.389264248198689e-4 ;F_t[2][2] = 0.999152490516521

                SHP = SHP3D(gaussw[0][ig],gaussw[1][ig],gaussw[2][ig],XYZ)
                F_tau = calc_F(xyz,SHP)

                F_tau[0][0] = 0.999966994751800 ; F_tau[0][1] = 6.616431928285574e-4 ; F_tau[0][2] = -4.680335771167036e-4
                F_tau[1][0] = 4.526431224534634e-7 ;  F_tau[1][1] = 1.00154297610866  ;  F_tau[1][2] = 2.184648048205991e-7
                F_tau[2][0] = -3.488804582905496e-6 ; F_tau[2][1] = 5.393786867065042e-4 ;  F_tau[2][2] = 0.999151405380578

########### Plastic Defomation Gradient at time = t should be here ######
##########  ??????????????????????????????
                Fp = [[0. for ii in range(3)] for jj in range(3)]
                Fp[0][0] = Fp[1][1] = Fp[2][2] = 1.0

                Fp[0][0] =  1.00001548434869 ;  Fp[0][1] = 3.336413447240032e-5 ; Fp[0][2] = 4.718401146770838e-5
                Fp[1][0] =  4.647041323286691e-5; Fp[1][1] = 1.00010012982441 ;  Fp[1][2] = 1.416049556808719e-4
                Fp[2][0] = -3.794103373053844e-5 ; Fp[2][1] = -8.175156580585659e-5 ; Fp[2][2] = 0.999884385826894
                
#------------------- Calculate matrix [A] 
                A_mat = calc_A(F_tau,Fp)
                
#------------------  Calculate matrix [B_alpha]
                B_alpha = calc_B(n_slip,A_mat,schmid)

#----------------- Calculate the trial elastic stress S_trial 
                S_trial = stress_trial(C_mat,A_mat)

#---------------- Calculate the matrix C_alpha for each slip system 
                C_alpha = calc_C_alpha(n_slip,C_mat,B_alpha)

######################## Enter the User Material ################### 
            
                res = [0.0 for ir in range(n_slip)]
                
                res[0] = 170846276.431964     
                res[1] = 170775924.340109     
                res[2] = 170846283.783709     
                res[3] = 170846283.783709     
                res[4] = 170846283.783709     
                res[5] = 170846283.783709     
                res[6] = 170846283.783709     
                res[7] = 170846283.783709     
                res[8] = 170846283.783709     
                res[9] = 170846283.783709     
                res[10] = 170846283.783709     
                res[11] = 170846283.783709     

                

                U = UMAT(S_trial,C_mat,Fp,C_alpha,schmid,F_t,F_tau,res)
#                U1 = U.itr()
                U2 = U.polardecomp()
                print U2[0]
                
    
#????????????????? End of UMAT 


#--------------- Updated Plastic Deforamtion Gradient at time tau 

#--------------- Updated Cauchy Stress 


#--------------- Updated Elastic and Plastic Strain

#                print "Element =",ielem.index+1,"Gauss point =",ig+1
#                print "Trial Stress"
#                for jslip  in range(3):
#                    print S_trial[jslip]

 #               for kslip in range(n_slip):
 #                   print "Slip System =",kslip+1
 #                   for islip  in range(3):
 #                       print C_alpha[kslip][islip]

            
        return  


m = Mesh()
m.solve()
#print m.solve()





























#########################   strain calculation ##################
def strain():
    icase = input("icase = ")
    if icase ==1:
        print "                  Stretch in X direction and Shrink in Y,Z"
    elif icase == 2:
        print "                          Rotation around Z axis"
    elif icase == 3:
        print "Stretch in X direction and Shrink in Y,Z then Rotation around Z axis"

        CORD = [[0. for ii in range(8)] for jj in range(3)]
        if icase == 1:
# x = 2X , y = 0.25Y , z = 0.25Z
            CORD[0] = [0.0,0.0,0.0,0.0,4.0,4.0,4.0,4.0]
            CORD[1] = [0.0,0.25,0.25,0.0,0.0,0.25,0.25,0.0]
            CORD[2] = [0.0,0.0,0.25,0.25,0.0,0.0,0.25,0.25]
        elif icase == 2:
# x = 0.8X-0.6Y , y = 0.6X+0.8Y , z = Z
            CORD[0] = [0.0,-0.6,-0.6,0.0,1.6,1.0,1.0,1.6]
            CORD[1] = [0.0,0.8,0.8,0.0,1.2,2.0,2.0,1.2]
            CORD[2] = [0.0,0.0,1.0,1.0,0.0,0.0,1.0,1.0]
        elif icase == 3:
# (x = 2X , y = 0.25Y , z = 0.25Z)*(x = 0.8X-0.6Y , y = 0.6X+0.8Y , z = Z)
            CORD[0] = [0.0,-0.15,-0.15,0.0,3.2,3.05,3.05,3.2]
            CORD[1] = [0.0,0.2,0.2,0.0,2.4,2.6,2.6,2.4]
            CORD[2] = [0.0,0.0,0.25,0.25,0.0,0.0,0.25,0.25]


#FG = [[[0. for ii in range(3)] for jj in range(3)] for kk in range(8)]

#    FGT=[]

#    gaussw = pgauss()


#    for ig in range(8):
#        XLT = [n1,n2,n3,n4,n5,n6,n7,n8]
#        XL = zip(*XLT)

#        SHP = SHP3D(gaussw[0][ig],gaussw[1][ig],gaussw[2][ig],XL)
#        F = CALC_F(CORD,SHP)
#        FGT.append(F)



    FAVE = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            for k in range(8):
                FAVE[i][j] += FGT[k][i][j]
                FAVE[i][j] = (1.0/8.0)*FAVE[i][j]


    FAVET = zip(*FAVE)
    Iden = [[0. for ii in range(3)] for jj in range(3)]
    Iden[0][0] = 1.0 ; Iden[1][1] = 1.0 ; Iden[2][2] = 1.0

    E = [[0. for ii in range(3)] for jj in range(3)]

# Green-Lagrange strain
    E = matmul(3,FAVET,FAVE)

    for i in range(3):
        for j in range(3):
            E[i][j] = 0.5*(E[i][j]-Iden[i][j])

    epsilonx = E[0][0] ; epsilonxy = E[0][1] ; epsilonxz = E[0][2]
    epsilony = E[1][1] ; epsilonyz = E[1][2] ; epsilonz = E[2][2]

    print "Green-Langrange strian components:"
    print "                                  Ex = %f" % epsilonx
    print "                                  Ey = %f" % epsilony
    print "                                  Ez = %f" % epsilonz
    print "                                  Exy = %f" % epsilonxy
    print "                                  Exz = %f" % epsilonxz
    print "                                  Eyz = %f" % epsilonyz


# Euler-Almansi strain
    e = [[0. for ii in range(3)] for jj in range(3)]
    e = matmul(3,FAVE,FAVET)

    eu = smallmatrix.SmallMatrix(3,3)
    Ieu = smallmatrix.SmallMatrix(3,3)
    eu.clear()
    Ieu.clear()
    for i in range(3):
        for j in range(3):
            Ieu[i,j] = 0.0
            Ieu[0,0]=Ieu[1,1]=Ieu[2,2]=1.0


    for i in range(3):
        for j in range(3):
            eu[i,j] = e[i][j]

    eu.solve(Ieu)

    for i in range(3):
        for j in range(3):
            e[i][j] = Ieu[i,j]


    for i in range(3):
        for j in range(3):
            e[i][j] = 0.5*(Iden[i][j]-e[i][j])

    epsilonx = e[0][0] ; epsilonxy = e[0][1] ; epsilonxz = e[0][2]
    epsilony = e[1][1] ; epsilonyz = e[1][2] ; epsilonz = e[2][2]

    print "Euler-Almansi strian components:"
    print "                                  ex = %f" % epsilonx
    print "                                  ey = %f" % epsilony
    print "                                  ez = %f" % epsilonz
    print "                                  exy = %f" % epsilonxy
    print "                                  exz = %f" % epsilonxz
    print "                                  eyz = %f" % epsilonyz


    du_dx = [[0. for ii in range(3)] for jj in range(3)]
    ep = [[0. for ii in range(3)] for jj in range(3)]

    for i in range(3):
        for j in range(3):
            du_dx[i][j] = FAVE[i][j]-Iden[i][j]
    for i in range(3):
        for j in range(3):
            ep[i][j] = 0.5*(du_dx[i][j]+du_dx[j][i])

    epsilonx = ep[0][0] ; epsilonxy = ep[0][1] ; epsilonxz = ep[0][2]
    epsilony = ep[1][1] ; epsilonyz = ep[1][2] ; epsilonz = ep[2][2]

    print "Small strian components:"
    print "                                  ex = %f" % epsilonx
    print "                                  ey = %f" % epsilony
    print "                                  ez = %f" % epsilonz
    print "                                  exy = %f" % epsilonxy
    print "                                  exz = %f" % epsilonxz
    print "                                  eyz = %f" % epsilonyz


####################### end of strain calculation ###########################


#for i in range(len(m.ellist)):
#    print m.ellist[i]
#print m.qrot
#print m.XYZ

#stress = [[0.0 for i in range(3)] for j in range(3)]
#stress[1][1] = 1.0
#tau_alpha = [0.0 for i in range(n_slip)]
#for k in range(n_slip):
#    for i in range(3):
#        for j in range(3):
#            tau_alpha[k] += stress[i][j]*s[k][i][j]

#for k in range(n_slip):
#    print ("%d %8.4f" % (k+1,tau_alpha[k]))

    
#for i in range(3):
#    for j in range(3):
#        for k in range(3):
#            for l in range(3):
#                print ("%d %d %d %d %12.4f" % (i+1,j+1,k+1,l+1,c_mat[i][j][k][l]))
