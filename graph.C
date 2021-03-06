#include <vector>

double sig_m = 1.;

double int_rho(double pos) {
	double sig = sig_m;
	double m = 1;
	double t = 5;
	double pi = TMath::Pi();
	TF2* rho = new TF2("rho", "(0.25/pi)*(sqrt(2./pi)*(1/[1]))*(1 + [3]*[3]/(sqrt([3]*[3] + x*x)*sqrt([3]*[3] + y*y)) + x*y/(sqrt([3]*[3] + x*x)*sqrt([3]*[3] + y*y)))*cos((x-y)*[0] - (sqrt([3]*[3] + x*x) - sqrt([3]*[3] + y*y))*[2])*TMath::Exp(-1.*(x*x + y*y)/([1]*[1])) ", -100*sig, 100*sig, -100*sig, 100*sig);
	rho->SetParameter(0, pos); //[0] == Position 
	rho->SetParameter(1, sig); //[1] == Standard deviation in momentum space
	rho->SetParameter(2, t);   //[2] == Time
	rho->SetParameter(3, m);   //[3] == Mass 
	return (rho->Integral(-5*sig, 5*sig, -5*sig, 5*sig,1.e-5)); 
}

double int_j(double pos) {
	double sig = sig_m;
	double m = 1;
	double t = 5;
	double pi = TMath::Pi();
	TF2* j = new TF2("j", "(0.25/pi)*(sqrt(2./pi)*(1/[1]))*(x/sqrt([3]*[3] + x*x) + y/sqrt([3]*[3] + y*y))*cos((x-y)*[0] - (sqrt([3]*[3] + x*x) - sqrt([3]*[3] + y*y))*[2])*TMath::Exp(-1.*(x*x + y*y)/([1]*[1])) ", -100*sig, 100*sig, -100*sig, 100*sig);
	j->SetParameter(0, pos);
	j->SetParameter(1, sig);
	j->SetParameter(2, t);
	j->SetParameter(3, m); 
	return (j->Integral(-5*sig, 5*sig, -5*sig, 5*sig, 1.e-5));

	//(sqrt(2)*sqrt(sqrt(2*pi)*[1]))/(([1]*sqrt([1]))
}


void graph() {
	//Defining rho and j as integrals
	double t = 5;
	double sig_0 = (1. / sig_m);
	//double sig = sig_0 * (1 + t * t / (4*sig_0 * sig_0));
	double sig = sig_0;
	double pi = TMath::Pi();
	double m = 1;
	TF1* g_rho = new TF1("g_rho", "int_rho(x)", -10*sig, 10*sig);
	TF1* g_j = new TF1("g_j", "int_j(x)", -10*sig, 10*sig);
	TF1* g_psi = new TF1("g_psi", "(1/sqrt(2*pi*[0]*[0]))*TMath::Exp(-0.5*x*x/([0]*[0]))", -20*sig, 20*sig);
	g_psi->SetParameter(0, sig);

	double rho_int = g_rho->Integral(-10 * sig, 10 * sig, 1.e-3);

//Attempt at faster 1D integration
/*	Int_t np = 1000;
	double *x_arr=new double[np];
       	double *w_arr=new double[np];
	g_rho->CalcGaussLegendreSamplingPoints(np,x_arr,w_arr,1e-15);
	double rho_int = g_rho->IntegralFast(np,x_arr,w_arr,-10*sig,10*sig);
*/
	//Normalisation check
	cout << "Integral of calculated density:" << rho_int << endl;
	cout << "Integral of classical density:" << g_psi->Integral(-10 * sig, 10 * sig, 1.e-3) << endl;

	//Renormalisation of j and rho (Ideally unnecessary)
	TF1* g_rho_fin = new TF1("g_rho_fin", "g_rho/[0]", -10 * sig, 10 * sig);
	TF1* g_j_fin = new TF1("g_j_fin", "g_j/[0]", -10 * sig, 10 * sig);
	g_rho_fin->SetParameter(0, rho_int);
	g_j_fin->SetParameter(0, rho_int);

	cout << "Renormalised integral of calculated density:" << g_rho_fin->Integral(-10 * sig, 10 * sig, 1.e-3) << endl;

	g_j->SetLineColor(1);
//	g_psi->SetLineColor(3);
	g_j_fin->SetLineColor(1);

	TCanvas* c1 = new TCanvas();
	c1->SetGrid();

	TMultiGraph* mg = new TMultiGraph();
	TGraph* g1 = new TGraph(g_rho_fin);
	TGraph* g2 = new TGraph(g_j_fin);
//	TGraph* g3 = new TGraph(g_psi);

	mg->Add(g1);
	mg->Add(g2);
//	mg->Add(g3);
	mg->Draw("apl");
	mg->GetXaxis()->SetTitle("x");
	mg->GetYaxis()->SetTitle("y");
	mg->GetHistogram()->SetTitle("Comparison of probability density with current for a Gaussian (sig=1,m=1,t=5)");

	TLegend* legend = new TLegend(0.78, 0.695, 0.98, 0.775);
	legend->AddEntry(g_rho_fin, "Probability density", "l");
	legend->AddEntry(g_j_fin, "Probability current", "l");
//	legend->AddEntry(g_psi, "Non-relativistic density", "l");
	legend->Draw("Same");
}

//This bit of code is for plotting the integrands and integral at a particular x (Needs mod for t)
/*void graph() {
	double pi = TMath::Pi();
	double pos = 1;
	double sig = 1;
	double m = 1;
	double t = 10;
	TF2* rho = new TF2("rho", "(1 + [3]*[3]/(sqrt([3]*[3] + x*x)*sqrt([3]*[3] + y*y)) + x*y/(sqrt([3]*[3] + x*x)*sqrt([3]*[3] + y*y)))*cos((x-y)*[0] - (sqrt([3]*[3] + x*x) - sqrt([3]*[3] + y*y))*[2])*TMath::Exp(-1.*(x*x + y*y)/([1]*[1])) ", -5*sig, 5*sig, -5*sig, 5*sig);
	TF2* j = new TF2("j", "(x/sqrt([3]*[3] + x*x) + y/sqrt([3]*[3] + y*y))*cos((x-y)*[0] - (sqrt([3]*[3] + x*x) - sqrt([3]*[3] + y*y))*[2])*TMath::Exp(-1.*(x*x + y*y)/([1]*[1])) ", -5*sig, 5*sig, -5*sig, 5*sig);
	rho->SetParameter(0, pos); //[0] == Position 
	rho->SetParameter(1, sig); //[1] == Standard deviation in momentum space
	rho->SetParameter(2, t);   //[2] == Time
	rho->SetParameter(3, m);   //[3] == Mass 
	j->SetParameter(0, pos);
	j->SetParameter(1, sig);
	j->SetParameter(2, t);
	j->SetParameter(3, m);
	double integral_rho =  0.5*(rho->Integral(-5 * sig, 5* sig, -5 * sig, 5 * sig));
	double integral_j = 0.5*(j->Integral(-5 * sig, 5* sig, -5 * sig, 5 * sig));


	TCanvas* c1 = new TCanvas();
	rho->Draw("colz");
	TCanvas* c2 = new TCanvas();
	j->Draw("colz");
	cout << "Probability density:" << integral_rho << endl;
	cout << "Probability current:" << integral_j << endl;
}*/
