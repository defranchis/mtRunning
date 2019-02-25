import ROOT
from ROOT import *
import os
import rundec


# input values (global variables)

xsec_1 = 255.37
err_xsec_1_up = 4.845
err_xsec_1_down = 4.497
extr_1_up = [ -0.152 , 2.666 , -0.117 , -0.025 , 2.134 , 2.326 , ]
extr_1_down = [ 0.679 , -0.183 , 0.449 , 0.239 , -1.428 , -3.133 , ]
err_xsec_toys_1 = 3.918

xsec_2 = 307.11
err_xsec_2_up = 4.864
err_xsec_2_down = 4.629
extr_2_up = [ -0.158 , 1.779 , -0.328 , -0.12 , 1.339 , -1.786 , ]
extr_2_down = [ 0.433 , -0.11 , 0.726 , 0.021 , -0.901 , 1.227 , ]
err_xsec_toys_2 = 5.134

xsec_3 = 180.86
err_xsec_3_up = 4.474
err_xsec_3_down = 4.234
extr_3_up = [ -0.21 , 1.605 , 0.648 , 0.134 , 1.019 , -2.629 , ]
extr_3_down = [ 0.046 , -0.042 , -0.114 , -0.079 , -0.844 , 1.804 , ]
err_xsec_toys_3 = 2.808

corr_1_2 = 0.7
corr_3_2 = 0.58

n_mttbins = 4

mass_min = 152.0
mass_max = 172.0
mass_fine_min = 152.0
mass_fine_max = 172.0

mu_1 = 383.95
mu_2 = 476.19
mu_3 = 644.27

MZ = 91.1876
asMZ = 0.11905
err_as = 0.0011
nflav = 5
nloops = 2
   
# end global variables


################################

# "setMasses" initializes which masses should be considered

################################

def setMasses():

    global mass_v
    mass_v = []
    i_mass = mass_max
    while i_mass >= mass_min :
        mass_v.append(i_mass)
        if i_mass <= mass_fine_max and i_mass > mass_fine_min+.1 : i_mass-= 0.2
        else : i_mass-= 1.0
    return


################################

# "mtmt2mtmu" converts mt(mt) to mt(mu) for a given scale
# as(MZ), MZ, n. of flavours and n. of loops should be set above

################################


def mtmt2mtmu(mt, mu):

    crd = rundec.CRunDec()

    asmt = crd.AlphasExact(asMZ, MZ, mt, nflav, nloops)
    asmu = crd.AlphasExact(asMZ, MZ, mu, nflav, nloops)
    
    mtmu = crd.mMS2mMS(mt, asmt, asmu, nflav, nloops)
    
    return mtmu



################################

# "mtmu2mtmu" converts mt(mu1) to mt(mu2) for a given scale
# as(MZ), MZ, n. of flavours and n. of loops should be set above

################################


def mtmu2mtmu (mtmu1, mu1, mu2, var_as):

    alphaS = 0
    if var_as == 'nominal' : alphaS = asMZ
    elif var_as == 'up' : alphaS = asMZ+err_as
    elif var_as == 'down' : alphaS = asMZ-err_as
    else :
        print 'ERROR!'
        return 0

    crd = rundec.CRunDec()
    
    asmu1 = crd.AlphasExact(alphaS, MZ, mu1, nflav, nloops)
    asmu2 = crd.AlphasExact(alphaS, MZ, mu2, nflav, nloops)
    
    mtmu2 = crd.mMS2mMS(mtmu1, asmu1, asmu2, nflav, nloops)
    
    return mtmu2


################################

# "formInputFileName" provides the correct name of the .dat input file from where
# the calculated cross section is read out

################################


def formInputFileName ( renscale, facscale, topmass, pdfmember ):

    infileName='tt_tot_tot_ABMP16_'
    infileName+=str(pdfmember)+'_'
    if pdfmember<10 : infileName+='_'
    
    if renscale>=1000 : infileName+='_'
    if renscale>=100 : infileName+=str(renscale)+'_'
    else : 
        renstr = TString(str(renscale))
        renstr.ReplaceAll('.','_.')
        infileName+=str(renstr)+'_'

    if facscale>=1000 : infileName+='_'
    if facscale>=100 : infileName+=str(facscale)+'_'
    else : 
        facstr = TString(str(facscale))
        facstr.ReplaceAll('.','_.')
        infileName+=str(facstr)+'_'

    infileName+=str(topmass)+'_MSbar.dat'

    #because of some problem with MCFM outputs...
    if TString(infileName).Contains('.2') :
        infileName=str(TString(infileName).ReplaceAll('.2','.1'))
    if TString(infileName).Contains('.6') :
        infileName=str(TString(infileName).ReplaceAll('.6','.5'))

        
    return infileName


################################

# "readCalculatedXsec" reads the calculated cross section in a given mtt bin
# for any value of mu_r, mu_f, mt and PDF member
# calls the function "formInputFileName"

################################


def readCalculatedXsec (renscale, facscale, topmass, pdfmember, mttbin):

    fileName = formInputFileName ( renscale, facscale, topmass, pdfmember )
    if not os.path.isfile('out_hists/'+fileName):
        print 'WARNING: missing file', fileName
        return 0
    infile = open('out_hists/'+fileName,'r')
    lines = infile.read().splitlines()
    fillList = False
    bincenters = []
    xsecs = []
    for line in lines: 
        if 'Cross-section' in line:
            xsec_tot = line.split(':')
            xsec_tot = xsec_tot[1].split('+/-')
            xsec_tot = float(xsec_tot[0])/1000.
        if 'HIST =   2' in line: fillList = True
        if fillList and not 'm34' in line and not 'HIST' in line:
            bin_center, bin_xsec, bin_err = line.split()
            bin_center = float(bin_center)
            bin_xsec = float(bin_xsec)
            if bin_xsec == 0 and bin_center > 300 : fillList = False
            bincenters.append(bin_center)
            xsecs.append(bin_xsec)

    bin_width=bincenters[1]-bincenters[0]

    if mttbin < 3:
        xsec = xsecs[mttbin]/1000.*bin_width
    elif n_mttbins == 3:
        # xsec=xsec_tot-(xsecs[1]+xsecs[2])/1000.*bin_width
        xsec=0
        for i in range(mttbin,mttbin+11):
            xsec += xsecs[i]/1000.*bin_width
    elif mttbin==3: # n_mttbins = 4
        xsec = (xsecs[3]+xsecs[4])/1000.*bin_width
    else:
        xsec=0
        for i in range(mttbin+1,mttbin+16):
            xsec += xsecs[i]/1000.*bin_width
            
        
    return xsec 


################################

# "getMassAndError" provides the extracted mass and its error in a given mtt bin
# and for a choice of mu_r, mu_f, pdf member and extrapolation uncertainty

################################

def getMassAndError(mttbin, murscale, mufscale, pdfmember, extrapol): 

    graph = ROOT.TGraph()
    i=-1
    chi2_v = []
    mass_v_sel = []
    for mass in mass_v:
        i+=1
        mur = mass
        muf = mass
        if (murscale == 'up')   : mur*=2
        if (mufscale == 'up')   : muf*=2
        if (murscale == 'down') : mur*=.5
        if (mufscale == 'down') : muf*=.5
        xsec_th = readCalculatedXsec(mur, muf, mass, pdfmember, mttbin) 
        if xsec_th == 0 : continue # torm
        xsec_exp = xsec_err = 0.
        if mttbin == 1 : 
            xsec_exp = xsec_1 
            xsec_err = ((xsec_1/100.*(err_xsec_1_up+err_xsec_1_down)/2.)**2 + err_xsec_toys_1**2)**.5
            if extrapol != 0:
                if extrapol > 0 :
                    xsec_exp *= 1 + extr_1_up[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_1_up[abs(extrapol)-1]/100.
                else :
                    xsec_exp *= 1 + extr_1_down[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_1_down[abs(extrapol)-1]/100.
        elif mttbin == 2 : 
            xsec_exp = xsec_2 
            xsec_err = ((xsec_2/100.*(err_xsec_2_up+err_xsec_2_down)/2.)**2 + err_xsec_toys_2**2)**.5
            if extrapol != 0:
                if extrapol > 0 :
                    xsec_exp *= 1 + extr_2_up[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_2_up[abs(extrapol)-1]/100.
                else :
                    xsec_exp *= 1 + extr_2_down[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_2_down[abs(extrapol)-1]/100.
        else : # mttbin = 3
            xsec_exp = xsec_3 
            xsec_err = ((xsec_3/100.*(err_xsec_3_up+err_xsec_3_down)/2.)**2 + err_xsec_toys_3**2)**.5
            if extrapol != 0:
                if extrapol > 0 :
                    xsec_exp *= 1 + extr_3_up[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_3_up[abs(extrapol)-1]/100.
                else :
                    xsec_exp *= 1 + extr_3_down[abs(extrapol)-1]/100.
                    xsec_err *= 1 + extr_3_down[abs(extrapol)-1]/100.


        # chi2 = abs(xsec_th-xsec_exp)/xsec_err
        chi2 = (xsec_th-xsec_exp)/xsec_err
        fact_A = 0
        if mttbin == 1 : fact_A = (err_xsec_1_up-err_xsec_1_down)/(err_xsec_1_up+err_xsec_1_down)
        elif mttbin == 2 : fact_A = (err_xsec_2_up-err_xsec_2_down)/(err_xsec_2_up+err_xsec_2_down)
        else : fact_A = (err_xsec_3_up-err_xsec_3_down)/(err_xsec_3_up+err_xsec_3_down)
        chi2 *= (1-2*fact_A*chi2+5*fact_A*fact_A*chi2*chi2)**.5
        # graph.SetPoint(i,mass,chi2)
        if mass != mass_v[0]:
            # if chi2<chi2_v[len(chi2_v)-1] : continue
            chi2_slope = (chi2-chi2_v[len(chi2_v)-1])/(mass-mass_v_sel[len(mass_v_sel)-1])
            slope_low = 0
            slope_high = 0
            if mttbin == 1 :
                slope_low = -6.0
                slope_high = -0.7
            elif mttbin == 2 :
                slope_low = -1.0
                slope_high = 0.5
            else:
                slope_low = -1.0
                slope_high = 0.05

            if chi2_slope < slope_low or chi2_slope > slope_high :
                if mass!= mass_v[1]: continue
                else:
                    if chi2 < chi2_v[len(chi2_v)-1] : continue

        mass_v_sel.append(mass)
        chi2_v.append(chi2)
            
    for i in range(0,len(chi2_v)):
        graph.SetPoint(i,mass_v_sel[i],chi2_v[i])

    funct = TF1('funct','pol4(0)',mass_min,mass_max)
    funct.SetParameter(0,1080)
    funct.SetParameter(1,-12)
    funct.SetParameter(2,0.033)
    # graph.Fit(funct,'q','',163,166)
    graph.Fit(funct,'q','',mass_min+0.1,mass_max-0.1)

    line_up = ROOT.TLine(mass_min,1.,mass_max,1.)
    line_down = ROOT.TLine(mass_min,-1.,mass_max,-1.)
    line_central = ROOT.TLine(mass_min,0.,mass_max,0.)
    line_up.SetLineColor(kGreen+2)
    line_down.SetLineColor(kGreen+2)
    line_central.SetLineColor(kRed)

    outdir = 'plots_chi2'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c = ROOT.TCanvas()
    c.SetGrid()
    graph.Draw('ap')
    graph.SetMarkerStyle(7)
    graph.SetTitle('signed #chi^{2} vs top mass in MCFM prediction; m_{t}(m_{t}) [GeV]; data-theory #chi^{2}')
    line_up.Draw('same')
    line_down.Draw('same')
    line_central.Draw('same')
    c.Print(outdir+'/test_mtt'+str(mttbin)+'_mur_'+murscale+'_muf_'+mufscale+'_extrapol'+str(extrapol)+'_pdf'+str(pdfmember)+'.png')

    fitted_mass = funct.GetX(0)
    fitted_mass_up = funct.GetX(-1)
    fitted_mass_down = funct.GetX(1)

    #now evolve masses
    mu = 0
    if mttbin == 1 : mu = mu_1
    elif mttbin == 2 : mu = mu_2
    else : mu = mu_3 # mttbin = 3

    if murscale=='nominal' and mufscale=='nominal' and pdfmember==0 and extrapol==0:
        print 'mt(mt) bin', mttbin,'=', fitted_mass, fitted_mass_up-fitted_mass, fitted_mass - fitted_mass_down
    
    fitted_mass = mtmt2mtmu(fitted_mass, mu)
    fitted_mass_up = mtmt2mtmu(fitted_mass_up, mu)
    fitted_mass_down = mtmt2mtmu(fitted_mass_down, mu)

    if murscale=='nominal' and mufscale=='nominal' and pdfmember==0 and extrapol==0:
        print 'mt(mu) bin', mttbin,'=', fitted_mass, fitted_mass_up-fitted_mass, fitted_mass - fitted_mass_down
        print
        
    fitted_mass_err = (fitted_mass_up - fitted_mass_down)/2

    return [fitted_mass, fitted_mass_err]


################################

# "getRatios" calculates all the ratios and their experimental errors
# taking the correlations into account

################################

def getRatios(mass_1, mass_2, mass_3, err_1, err_2, err_3):

    ratio_1_2 = mass_1/mass_2
    err_ratio_1_2 = (err_1/mass_1)**2 + (err_2/mass_2)**2 - 2*corr_1_2*(err_1/mass_1)*(err_2/mass_2)
    err_ratio_1_2 = err_ratio_1_2**.5
    err_ratio_1_2*= ratio_1_2

    ratio_3_2 = mass_3/mass_2
    err_ratio_3_2 = (err_3/mass_3)**2 + (err_2/mass_2)**2 - 2*corr_3_2*(err_3/mass_3)*(err_2/mass_2)
    err_ratio_3_2 = err_ratio_3_2**.5
    err_ratio_3_2*= ratio_3_2

    return [ratio_1_2, ratio_3_2, err_ratio_1_2, err_ratio_3_2]
    

################################

# "getScaleUncertainties" calculates the scale uncertainties of the ratios

################################


def getScaleUncertainties (central_ratio_1_2, central_ratio_3_2):

    err_scale_up_1_2 = central_ratio_1_2
    err_scale_up_3_2 = central_ratio_3_2
    err_scale_down_1_2 = central_ratio_1_2
    err_scale_down_3_2 = central_ratio_3_2

    for renscale in 'nominal', 'up','down':
        for facscale in 'nominal', 'up','down':

            if renscale == 'up' and facscale == 'down': continue
            if facscale == 'up' and renscale == 'down': continue
            if facscale == 'nominal' and renscale == 'nominal': continue
            
            mass_and_err_1 = getMassAndError(1, renscale, facscale, 0 , 0)
            mass_and_err_2 = getMassAndError(2, renscale, facscale, 0 , 0)
            mass_and_err_3 = getMassAndError(3, renscale, facscale, 0 , 0)

            ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0],
                                        mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1])

            ratio_1_2 = ratios_and_errs[0]
            ratio_3_2 = ratios_and_errs[1]
    
            if ratio_1_2 > err_scale_up_1_2 : err_scale_up_1_2 = ratio_1_2
            if ratio_3_2 > err_scale_up_3_2 : err_scale_up_3_2 = ratio_3_2
            if ratio_1_2 < err_scale_down_1_2 : err_scale_down_1_2 = ratio_1_2
            if ratio_3_2 < err_scale_down_3_2 : err_scale_down_3_2 = ratio_3_2
    

    err_scale_up_1_2 -= central_ratio_1_2
    err_scale_up_3_2 -= central_ratio_3_2
    err_scale_down_1_2 = central_ratio_1_2 - err_scale_down_1_2
    err_scale_down_3_2 = central_ratio_3_2 - err_scale_down_3_2

    return [err_scale_up_1_2, err_scale_down_1_2, err_scale_up_3_2, err_scale_down_3_2]
    

################################

# "getPDFUncertainties" calculates the PDF uncertainties of the ratios

################################

def getPDFUncertainties (central_ratio_1_2, central_ratio_3_2):

    err_pdf_up_1_2 = 0
    err_pdf_up_3_2 = 0
    err_pdf_down_1_2 = 0
    err_pdf_down_3_2 = 0
    
    for pdf in range(1,30):
        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', pdf , 0)
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', pdf , 0)
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', pdf , 0)

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1])

        ratio_1_2 = ratios_and_errs[0]
        ratio_3_2 = ratios_and_errs[1]

        err_pdf_1_2 = (ratio_1_2-central_ratio_1_2)**2
        err_pdf_3_2 = (ratio_3_2-central_ratio_3_2)**2

        if ratio_1_2 > central_ratio_1_2 : err_pdf_up_1_2 += err_pdf_1_2
        else : err_pdf_down_1_2 += err_pdf_1_2

        if ratio_3_2 > central_ratio_3_2 : err_pdf_up_3_2 += err_pdf_3_2
        else : err_pdf_down_3_2 += err_pdf_3_2
        
    return [err_pdf_up_1_2**.5, err_pdf_down_1_2**.5, err_pdf_up_3_2**.5, err_pdf_down_3_2**.5]
        

################################

# "getExtrapolationUncertainties" calculates the extrapolation uncertainties of the ratios

################################

def getExtrapolationUncertainties (central_ratio_1_2, central_ratio_3_2):

    err_extr_up_1_2 = 0
    err_extr_up_3_2 = 0
    err_extr_down_1_2 = 0
    err_extr_down_3_2 = 0

    for extr in range(-len(extr_1_up),len(extr_1_up)+1) :
        if extr == 0 : continue

        mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , extr)
        mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , extr)
        mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , extr)

        ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0],
                                    mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1])

        ratio_1_2 = ratios_and_errs[0]
        ratio_3_2 = ratios_and_errs[1]

        err_extr_1_2 = (ratio_1_2-central_ratio_1_2)**2
        err_extr_3_2 = (ratio_3_2-central_ratio_3_2)**2

        if ratio_1_2 > central_ratio_1_2 : err_extr_up_1_2 += err_extr_1_2
        else : err_extr_down_1_2 += err_extr_1_2

        if ratio_3_2 > central_ratio_3_2 : err_extr_up_3_2 += err_extr_3_2
        else : err_extr_down_3_2 += err_extr_3_2


    return [err_extr_up_1_2**.5, err_extr_down_1_2**.5, err_extr_up_3_2**.5, err_extr_down_3_2**.5]


################################

# "makeTheoryPrediction" calculates the running of mt(mu) starting from m(mu_2)

################################


def makeTheoryPrediction(outfile, mass_2):

    scales = []
    r_up = []
    r_down = []
    out = open(outfile+'.txt','w')
        
    for scale in range(350,670+1):
        ratio = mtmu2mtmu(mass_2, mu_2, scale, 'nominal')/mass_2
        out.write(str(scale)+'\t'+str(ratio)+'\n')
        ratio_up = mtmu2mtmu(mass_2, mu_2, scale, 'up')/mass_2
        ratio_down = mtmu2mtmu(mass_2, mu_2, scale, 'down')/mass_2
        r_up.append(ratio_up)
        r_down.append(ratio_down)
        scales.append(scale)
        
    out.close()
        
    return [scales, r_up, r_down]



################################

# "makeAdditionalTheoryPrediction" calculates the running of mt(mu) starting from m(mt)

################################


def makeAdditionalTheoryPrediction (mtmt, err_mtmt, mtmu):
    r = []
    ru = []
    rd = []
    scales = []
    r.append(mtmt/mtmu)
    ru.append((mtmt+err_mtmt)/mtmu)
    rd.append((mtmt-err_mtmt)/mtmu)
    scales.append(mtmt)
    for scale in range(int(mtmt)+1,670+1):
        ratio = mtmt2mtmu(mtmt,scale)/mtmu
        ratio_up = mtmt2mtmu(mtmt+err_mtmt,scale)/mtmu
        ratio_down = mtmt2mtmu(mtmt-err_mtmt,scale)/mtmu
        r.append(ratio)
        ru.append(ratio_up)
        rd.append(ratio_down)
        scales.append(scale)
        
    return [r, ru, rd, scales]



################################

# "makePlots" produces the plots for the running of mt

################################


def makePlots (mass_2, ratio_12, ratio_32, err_12, err_32):

    graph = TGraphErrors(3)

    graph.SetPoint(0,mu_1,ratio_12)
    graph.SetPointError(0,0,err_12)

    graph.SetPoint(1,mu_2,1)
    graph.SetPointError(1,0,0)

    graph.SetPoint(2,mu_3,ratio_32)
    graph.SetPointError(2,0,err_32)

    theoryFileName = 'theory_prediction'
    l = makeTheoryPrediction(theoryFileName,mass_2)
    th = TGraph(theoryFileName+'.txt')
    th.SetLineColor(kRed);
    th.SetLineWidth(2);

    scales = l[0]
    r_up = l[1]
    r_down = l[2]
    
    th_band = TGraph (2*th.GetN())

    for i in range(0, th.GetN()):
        th_band.SetPoint(i,scales[i],max(r_up[i],r_down[i]))
        th_band.SetPoint(th.GetN()+i,scales[th.GetN()-i-1],min(r_up[th.GetN()-i-1],r_down[th.GetN()-i-1]));

    th_band.SetFillStyle(3001);
    th_band.SetFillColor(kRed);
    
    leg = TLegend(.15,.2,.6,.32)
    leg.AddEntry(graph,'MCFM @NLO from diff. #sigma_{t#bar{t}}','pe')
    leg.AddEntry(th,'RunDec @ 2 loops (5 flav.)','l')

    latexLabel1 = TLatex();
    latexLabel1.SetTextSize(0.04);
    latexLabel1.SetTextFont(42)
    latexLabel1.SetNDC();
    
    graph.GetXaxis().SetTitle('median m_{t#bar{t}} [GeV]')
    graph.GetYaxis().SetTitle('m_{t}(m_{t#bar{t}}) / m_{t}('+str(int(mu_2))+' GeV)')
    graph.SetTitle('running of m_{t}(#mu) as a function of #mu=m_{t#bar{t}}')

    c = TCanvas()
    graph.SetMarkerStyle(8)
    graph.Draw('ap')
    th.Draw("L same");
    leg.Draw('same')
    th_band.Draw('f same')
    latexLabel1.DrawLatex(0.59, 0.78, "ABMP16_5_nlo PDF set");

    outdir = 'plots_running'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c.SaveAs(outdir+'/running.png')
    c.SaveAs(outdir+'/running.pdf')
    c.SaveAs(outdir+'/running.root')

    mtmt = 164.02
    mtmt_err = 1.58

    mtmu = mtmt2mtmu (mtmt, mu_2)

    gr_add = TGraphErrors()

    gr_add.SetPoint(0,mtmt,mtmt/mtmu)
    gr_add.SetPointError(0,0,mtmt_err/mtmu)

    l = makeAdditionalTheoryPrediction(mtmt, mtmt_err, mtmu)
    r = l[0]
    r_up = l[1]
    r_down = l[2]
    scales = l[3]

    gr_band = TGraph (2*gr_add.GetN())

    for i in range(0, len(r)):
        gr_band.SetPoint(i,scales[i],r_up[i])
        gr_band.SetPoint(len(r)+i,scales[len(r)-i-1],r_down[len(r)-i-1]);

    gr_band.SetFillStyle(3002);
    gr_band.SetFillColor(kRed+1);

    gr_band.GetXaxis().SetTitle('#mu [GeV]')
    gr_band.GetYaxis().SetTitle('m_{t}(#mu) / m_{t}('+str(int(mu_2))+' GeV)')
    gr_band.SetTitle('running of m_{t}(#mu) as a function of #mu')

    gr_add.SetMarkerStyle(8)
    gr_add.SetMarkerColor(kBlue)
    gr_add.SetLineColor(kBlue)

    gr_band.GetYaxis().SetRangeUser(0.9,1.1)

    leg2 = leg.Clone()
    leg2.Clear()
    leg2.AddEntry(graph,'MCFM @NLO from diff. #sigma_{t#bar{t}}','pe')
    leg2.AddEntry(gr_add,'Hathor @NLO from incl. #sigma_{t#bar{t}}')
    leg2.AddEntry(gr_band,'RunDec @ 2 loops (5 flav.)','f')

    g = graph.Clone()
    g.RemovePoint(1)
    g1 = TGraph()
    g1.SetPoint(0,mu_2,1)
    g1.SetMarkerStyle(3)
    g1.SetMarkerSize(1.5)
    
    c.Clear()
    gr_band.Draw('af')
    gr_add.Draw('p same')
    g.Draw('p same')
    g1.Draw('p same')
    # th.Draw("L same");
    leg2.Draw('same')
    # th_band.Draw('f same')

    c.SaveAs(outdir+'/test_incl.png')
    c.SaveAs(outdir+'/test_incl.pdf')
    c.SaveAs(outdir+'/test_incl.root')

    return


################################

# main program

################################


def execute():

    setMasses()
    mass_and_err_1 = getMassAndError(1, 'nominal', 'nominal', 0 , 0)
    mass_and_err_2 = getMassAndError(2, 'nominal', 'nominal', 0 , 0)
    mass_and_err_3 = getMassAndError(3, 'nominal', 'nominal', 0 , 0)

    ratios_and_errs = getRatios(mass_and_err_1[0], mass_and_err_2[0], mass_and_err_3[0],
                                mass_and_err_1[1], mass_and_err_2[1], mass_and_err_3[1])

    ratio_1_2 = ratios_and_errs[0]
    ratio_3_2 = ratios_and_errs[1]
    err_ratio_1_2 = ratios_and_errs[2]
    err_ratio_3_2 = ratios_and_errs[3]


    pdf_errors = getPDFUncertainties(ratio_1_2, ratio_3_2)

    err_pdf_1_2_up = pdf_errors[0]
    err_pdf_1_2_down = pdf_errors[1]
    err_pdf_3_2_up = pdf_errors[2]
    err_pdf_3_2_down = pdf_errors[3]


    extr_errors = getExtrapolationUncertainties(ratio_1_2, ratio_3_2)

    err_extr_1_2_up = extr_errors[0]
    err_extr_1_2_down = extr_errors[1]
    err_extr_3_2_up = extr_errors[2]
    err_extr_3_2_down = extr_errors[3]

    err_1_2_up = (err_ratio_1_2**2 + err_pdf_1_2_up**2 + err_extr_1_2_up **2)**.5
    err_1_2_down = (err_ratio_1_2**2 + err_pdf_1_2_down**2 + err_extr_1_2_down **2)**.5

    err_3_2_up = (err_ratio_3_2**2 + err_pdf_3_2_up**2 + err_extr_3_2_up **2)**.5
    err_3_2_down = (err_ratio_3_2**2 + err_pdf_3_2_down**2 + err_extr_3_2_down **2)**.5

    print 'ratio_1_2 =', ratio_1_2, '+' , err_1_2_up, '-' , err_1_2_down 
    print 'ratio_3_2 =', ratio_3_2, '+' , err_3_2_up, '-' , err_3_2_down 
    print

    print 'ratio_1_2 =', round(ratio_1_2,3), '+' , round(err_1_2_up,3), '-' , round(err_1_2_down,3)
    print 'ratio_3_2 =', round(ratio_3_2,3), '+' , round(err_3_2_up,3), '-' , round(err_3_2_down,3) 
    print

    makePlots (mass_and_err_2[0], ratio_1_2, ratio_3_2, err_ratio_1_2, err_ratio_3_2)

    return


################################

# execute main program

################################


execute()









