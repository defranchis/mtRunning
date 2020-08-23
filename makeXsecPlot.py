import ROOT
from ROOT import TString, TFile, TGraph, TGraphAsymmErrors, TCanvas
import os

preliminary = False
cms = True

th_xsec_164 = [261.4, 304.4, 177.3, 46.6]
th_xsec_162 = [294.1, 314.8, 181.0, 47.3]
th_xsec_166 = [231.4, 294.7, 173.3, 46.2]

tot_164 = th_xsec_164[0] + th_xsec_164[1] + th_xsec_164[2] + th_xsec_164[3]
tot_162 = th_xsec_162[0] + th_xsec_162[1] + th_xsec_162[2] + th_xsec_162[3]
tot_166 = th_xsec_166[0] + th_xsec_166[1] + th_xsec_166[2] + th_xsec_166[3]


bins = [305,420,550,810,1995]
scales = [384.0, 476.2, 644.3, 1023.6]


def getXsec(mttbin):

    diff = [0.318451076954, 0.392863602995, 0.226057763342, 0.0626275567086]
    err = [0.00687182481433, 0.0075685508001, 0.0056250500562, 0.0030817548822]
    
    return [diff[mttbin-1], err[mttbin-1], err[mttbin-1]]


g_th_164 = TGraphAsymmErrors()
g_th_162 = TGraphAsymmErrors()
g_th_166 = TGraphAsymmErrors()
g_exp = TGraphAsymmErrors()

for bin in range(1,5):
    g_th_164.SetPoint(bin-1, scales[bin-1], th_xsec_164[bin-1]/tot_164)
    g_th_164.SetPointError(bin-1,scales[bin-1]-bins[bin-1], bins[bin]-scales[bin-1], 0, 0)
    g_th_162.SetPoint(bin-1, scales[bin-1], th_xsec_162[bin-1]/tot_162)
    g_th_162.SetPointError(bin-1,scales[bin-1]-bins[bin-1], bins[bin]-scales[bin-1], 0, 0)
    g_th_166.SetPoint(bin-1, scales[bin-1], th_xsec_166[bin-1]/tot_166)
    g_th_166.SetPointError(bin-1,scales[bin-1]-bins[bin-1], bins[bin]-scales[bin-1], 0, 0)

    xsec_err = getXsec(bin)
    g_exp.SetPoint(bin-1, scales[bin-1], xsec_err[0])
    g_exp.SetPointError(bin-1, 0, 0, xsec_err[2], xsec_err[1])


g_th_164.SetMarkerStyle(0)
g_th_164.SetMarkerColor(ROOT.kRed)
g_th_164.SetLineColor(ROOT.kRed)

g_th_162.SetMarkerStyle(0)
g_th_162.SetMarkerColor(ROOT.kBlue)
g_th_162.SetLineColor(ROOT.kBlue)
g_th_162.SetLineStyle(2)

g_th_166.SetMarkerStyle(0)
g_th_166.SetMarkerColor(ROOT.kGreen+2)
g_th_166.SetLineColor(ROOT.kGreen+2)
g_th_166.SetLineStyle(3)

g_exp.SetMarkerStyle(8)
g_exp.SetLineWidth(1)


latexLabel1 = ROOT.TLatex()
latexLabel1.SetTextSize(0.06)
latexLabel1.SetNDC()

latexLabel2 = ROOT.TLatex()
latexLabel2.SetTextSize(0.04)
latexLabel2.SetTextFont(42)
latexLabel2.SetNDC()

latexLabel3 = ROOT.TLatex()
latexLabel3.SetTextSize(0.04)
latexLabel3.SetTextFont(42)
latexLabel3.SetNDC()

leg1 = ROOT.TLegend(.40,.75,.89,.85)
leg1.SetBorderSize(0)
leg1.AddEntry(g_exp,'Data unfolded to parton level','pe')
leg2 = ROOT.TLegend(.52,.33,.89,.48)
leg2.SetBorderSize(0)
leg2.AddEntry(g_th_162,'m_{t}(m_{t}) = 162 GeV','pl')
leg2.AddEntry(g_th_164,'m_{t}(m_{t}) = 164 GeV','pl')
leg2.AddEntry(g_th_166,'m_{t}(m_{t}) = 166 GeV','pl')

c = TCanvas()
c.SetLeftMargin(0.14)
c.SetBottomMargin(0.14)
c.SetRightMargin(0.04)
c.SetTopMargin(0.08)
g_th_162.GetXaxis().SetTitle('m_{t#bar{t}} [GeV]')
g_th_162.GetYaxis().SetTitle('#frac{1}{#sigma_{t#bar{t}}} d#sigma_{t#bar{t}} / dm_{t#bar{t}} #Deltam_{t#bar{t}}')
g_th_162.GetYaxis().SetTitleSize(0.05)
g_th_162.GetYaxis().SetTitleOffset(1.2)
g_th_162.GetXaxis().SetTitleSize(0.06)
g_th_162.GetXaxis().SetTitleOffset(1)
g_th_162.GetYaxis().SetLabelSize(0.05)
g_th_162.GetXaxis().SetLabelSize(0.05)
g_th_162.GetXaxis().SetRangeUser(bins[0],bins[len(bins)-1])
g_th_162.GetXaxis().SetLimits(137,2162)
g_th_162.SetMaximum(.43)


g_th_162.Draw('ap')
g_th_162.Draw('psame')
g_exp.Draw('psame')
g_th_164.Draw('psame')
g_th_166.Draw('psame')
g_exp.Draw('psame')
leg1.Draw('same')
leg2.Draw('same')


if cms: latexLabel1.DrawLatex(0.16, 0.94, "CMS")
latexLabel2.DrawLatex(0.76, 0.94, "35.9 fb^{-1} (13 TeV)")
if preliminary: latexLabel3.DrawLatex(0.215, 0.92 , "#it{Preliminary}")

latexLabel3.DrawLatex(.523, 0.67, "NLO predictions in #bar{MS} scheme")
latexLabel3.DrawLatex(.523, 0.62, "#mu_{r} = #mu_{f} = m_{t}")
latexLabel3.DrawLatex(.523, 0.54, "ABMP16_5_nlo PDF set")

outdir = 'xsec_plots'

if not os.path.exists(outdir):
    os.makedirs(outdir)



c.SaveAs(outdir+'/diffXsec.png')
c.SaveAs(outdir+'/diffXsec.pdf')
c.SaveAs(outdir+'/diffXsec.root')

