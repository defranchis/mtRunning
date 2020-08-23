import variables as var
from variables import xsec_1, xsec_2, xsec_3, xsec_4

import os

import ROOT as rt
from ROOT import TRandom3, TH2D, TH1D, TCanvas

import numpy as np

ntoys = 10000

def throwToyCrossSections(r):

    err_1 = xsec_1*(var.err_xsec_1_up+var.err_xsec_1_down)/2./100.
    err_2 = xsec_2*(var.err_xsec_2_up+var.err_xsec_2_down)/2./100.
    err_3 = xsec_3*(var.err_xsec_3_up+var.err_xsec_3_down)/2./100.
    err_4 = xsec_4*(var.err_xsec_4_up+var.err_xsec_4_down)/2./100.

    V = [[err_1*err_1, var.corr_1_2*err_1*err_2, var.corr_1_3*err_1*err_3, var.corr_4_1*err_1*err_4],
         [var.corr_1_2*err_1*err_2, err_2*err_2, var.corr_3_2*err_2*err_3, var.corr_4_2*err_2*err_4],
         [var.corr_1_3*err_1*err_3, var.corr_3_2*err_2*err_3, err_3*err_3, var.corr_4_3*err_3*err_4],
         [var.corr_4_1*err_1*err_4, var.corr_4_2*err_2*err_4, var.corr_4_3*err_3*err_4, err_4*err_4]]

    L = np.linalg.cholesky(V) # Cholesky decomposition 

    mu = [xsec_1, xsec_2, xsec_3, xsec_4]
    z = [r.Gaus(0,1),r.Gaus(0,1),r.Gaus(0,1),r.Gaus(0,1)]

    y = L.dot(z)
    y += mu

    return y


def execute():

    r=TRandom3()
    
    min_1 = .28
    max_1 = .36
    min_2 = .35
    max_2 = .44
    min_3 = .19
    max_3 = .26
    min_4 = .045
    max_4 = .08
    
    c_12 = TH2D('c_12','c_12',100,min_1,max_1,100,min_2,max_2)
    c_13 = TH2D('c_13','c_13',100,min_1,max_1,100,min_3,max_3)
    c_14 = TH2D('c_14','c_14',100,min_1,max_1,100,min_4,max_4)
    c_23 = TH2D('c_23','c_23',100,min_2,max_2,100,min_3,max_3)
    c_24 = TH2D('c_24','c_24',100,min_2,max_2,100,min_4,max_4)
    c_34 = TH2D('c_34','c_34',100,min_3,max_3,100,min_4,max_4)

    incl = TH1D('incl','incl',100,600,1000)
    
    for i in range(1,ntoys+1):
        if i%10000==0 or i==1: print 'toy n.', i
        res = throwToyCrossSections(r)

        tot = res[0] + res[1] + res[2] + res[3]
        diff_1 = res[0]/tot
        diff_2 = res[1]/tot
        diff_3 = res[2]/tot
        diff_4 = res[3]/tot
        c_12.Fill(diff_1, diff_2)
        c_13.Fill(diff_1, diff_3)
        c_14.Fill(diff_1, diff_4)
        c_23.Fill(diff_2, diff_3)
        c_24.Fill(diff_2, diff_4)
        c_34.Fill(diff_3, diff_4)
        incl.Fill(tot)
        # print tot, diff_1, diff_2, diff_3, diff_4

    outdir = 'corr_plots/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)


    c = TCanvas()
    c_12.SetTitle('; bin 1; bin 2')
    c_12.SetTitle('corr = '+str(round(c_12.GetCorrelationFactor(),2)))
    c_12.DrawNormalized('colz')
    c.SaveAs(outdir+'c_12.png','png')
    c.Clear()
    
    c_13.SetTitle('; bin 1; bin 3')
    c_13.SetTitle('corr = '+str(round(c_13.GetCorrelationFactor(),2)))
    c_13.DrawNormalized('colz')
    c.SaveAs(outdir+'c_13.png','png')
    c.Clear()

    c_14.SetTitle('; bin 1; bin 4')
    c_14.SetTitle('corr = '+str(round(c_14.GetCorrelationFactor(),2)))
    c_14.DrawNormalized('colz')
    c.SaveAs(outdir+'c_14.png','png')
    c.Clear()
    
    c_23.SetTitle('; bin 2; bin 3')
    c_23.SetTitle('corr = '+str(round(c_23.GetCorrelationFactor(),2)))
    c_23.DrawNormalized('colz')
    c.SaveAs(outdir+'c_23.png','png')
    c.Clear()
    
    c_24.SetTitle('; bin 2; bin 4')
    c_24.SetTitle('corr = '+str(round(c_24.GetCorrelationFactor(),2)))
    c_24.DrawNormalized('colz')
    c.SaveAs(outdir+'c_24.png','png')
    c.Clear()

    c_34.SetTitle('; bin 3; bin 4')
    c_34.SetTitle('corr = '+str(round(c_34.GetCorrelationFactor(),2)))
    c_34.DrawNormalized('colz')
    c.SaveAs(outdir+'c_34.png','png')
    c.Clear()


    incl.DrawNormalized()
    c.SaveAs(outdir+'incl.png','png')


    diff_1 = c_12.ProjectionX().GetMean()
    err_1 = c_12.ProjectionX().GetRMS()
    diff_2 = c_12.ProjectionY().GetMean()
    err_2 = c_12.ProjectionY().GetRMS()
    diff_3 = c_34.ProjectionX().GetMean()
    err_3 = c_34.ProjectionX().GetRMS()
    diff_4 = c_34.ProjectionY().GetMean()
    err_4 = c_34.ProjectionY().GetRMS()

    c12 = c_12.GetCovariance()
    c13 = c_13.GetCovariance()
    c14 = c_14.GetCovariance()
    c23 = c_23.GetCovariance()
    c24 = c_24.GetCovariance()
    c34 = c_34.GetCovariance()
    

    V = [[err_1**2, c12, c13, c14],
         [c12, err_2**2, c23, c24],
         [c13, c23, err_3**2, c34],
         [c14, c24, c34, err_4**2],]
    
    x = [diff_1, diff_2, diff_3, diff_4]

    print diff_1, diff_2, diff_3, diff_4
    print err_1, err_2, err_3, err_4
        
    
    return


execute()
