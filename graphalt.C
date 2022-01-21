#include <vector>

double int_rho(double pos) {
	double sig = 1.;
	double m = 1;
	double t = 0;
	TF2* rho = new TF2("rho", "(1 + [3]*[3]/(sqrt([3]*[3] + x*x)*sqrt([3]*[3] + y*y)) + x*y/(sqrt([3]*[3] + x*x)*sqrt([3]*[3] + y*y)))*cos((x-y)*[0] - (sqrt([3]*[3] + x*x) - sqrt([3]*[3] + y*y))*[2])*TMath::Exp(-1.*(x*x + y*y)/([1]*[1])) ", -1, 1, -1, 1);
	rho->SetParameter(0, pos); //[0] == Position 
	rho->SetParameter(1, sig); //[1] == Standard deviation in momentum space
	rho->SetParameter(2, t);   //[2] == Time
	rho->SetParameter(3, m);   //[3] == Mass 
	TF2* rho_n = new TF2("rho_n", "(1 + ([3]*[3]*x*y)/(sqrt(x*x*[3]*[3] + 1)*sqrt(y*y*[3]*[3] + 1)) + 1/(sqrt(x*x*[3]*[3] + 1)*sqrt(y*y*[3]*[3] + 1)))*cos([0]*(y-x)/(x*y) - (sqrt(x*x*[3]*[3] + 1)/x - sqrt(y*y*[3]*[3] + 1)/y)*[2])*TMath::Exp(-1.*(x*x + y*y)/([1]*[1]*x*x*y*y))*(1/(x*x*y*y)) ", -1, 1, -1, 1);
	rho_n->SetParameter(0, pos); //[0] == Position 
	rho_n->SetParameter(1, sig); //[1] == Standard deviation in momentum space
	rho_n->SetParameter(2, t);   //[2] == Time
	rho_n->SetParameter(3, m);   //[3] == Mass 
	return 0.5 * (rho->Integral(0, 1, 0, 1, 1.e-6) + rho_n->Integral(0,1,0,1,1.e-6) + rho->Integral(-1, 0, -1, 0, 1.e-6) + rho_n->Integral(-1, 0, -1, 0, 1.e-6));
}

double int_j(double pos) {
	double sig = 1.;
	double m = 1;
	double t = 0;
	TF2* j = new TF2("j", "(x/sqrt([3]*[3] + x*x) + y/sqrt([3]*[3] + y*y))*cos((x-y)*[0] - (sqrt([3]*[3] + x*x) - sqrt([3]*[3] + y*y))*[2])*TMath::Exp(-1.*(x*x + y*y)/([1]*[1])) ", -1, 1, -1, 1);
	j->SetParameter(0, pos);
	j->SetParameter(1, sig);
	j->SetParameter(2, t);
	j->SetParameter(3, m);
	TF2* j_n = new TF2("j_n", "(1/sqrt(x*x*[3]*[3] + 1) + 1/sqrt(y*y*[3]*[3] + 1))*cos([0]*(y-x)/(x*y) - (sqrt(x*x*[3]*[3] + 1)/x - sqrt(y*y*[3]*[3] + 1)/y)*[2])*TMath::Exp(-1.*(x*x + y*y)/([1]*[1]*x*x*y*y))*(1/(x*x*y*y))", -1, 1, -1, 1);
	j_n->SetParameter(0, pos);
	j_n->SetParameter(1, sig);
	j_n->SetParameter(2, t);
	j_n->SetParameter(3, m);
	return 0.5 * (j->Integral(0, 1, 0, 1, 1.e-6) + j_n->Integral(0, 1, 0, 1, 1.e-6) + j->Integral(-1, 0, -1, 0, 1.e-6) + j_n->Integral(-1, 0, -1, 0, 1.e-6));

	//(sqrt(2)*sqrt(sqrt(2*pi)*[1]))/(([1]*sqrt([1]))
}



void graphalt() {
	//Defining rho and j as integrals
	double t = 0;
	double sig_0 = (1. / 1.);
	//double sig = sig_0 * (1 + t * t / (4*sig_0 * sig_0));
	double sig = sig_0;
	double pi = TMath::Pi();
	double m = 1;
	TF1* g_rho = new TF1("g_rho", "int_rho(x)/(2*pi)", -10 * sig, 10 * sig);
	TF1* g_j = new TF1("g_j", "(int_j(x)/(2*pi))", -10 * sig, 10 * sig);
	TF1* g_psi = new TF1("g_psi", "(1/sqrt(2*pi*[0]*[0]))*TMath::Exp(-0.5*x*x/([0]*[0]))", -10 * sig, 10 * sig);
	g_psi->SetParameter(0, sig);

	double rho_int = g_rho->Integral(-10 * sig, 10 * sig, 1.e-3);

	//Normalisation check
	cout << "Integral of calculated density:" << rho_int << endl;
	cout << "Integral of classical density:" << g_psi->Integral(-10 * sig, 10 * sig, 1.e-3) << endl;

	//Renormalisation of j and rho (Ideally unnecessary, but clearly still needed)
	TF1* g_rho_fin = new TF1("g_rho_fin", "g_rho/[0]", -10 * sig, 10 * sig);
	TF1* g_j_fin = new TF1("g_j_fin", "g_j/[0]", -10 * sig, 10 * sig);
	g_rho_fin->SetParameter(0, rho_int);
	g_j_fin->SetParameter(0, rho_int);

	cout << "Renormalised integral of calculated density:" << g_rho_fin->Integral(-10 * sig, 10 * sig, 1.e-3) << endl;

	g_j->SetLineColor(1);
	g_psi->SetLineColor(3);
	g_j_fin->SetLineColor(1);

	TCanvas* c1 = new TCanvas();
	c1->SetGrid();

	TMultiGraph* mg = new TMultiGraph();
	TGraph* g1 = new TGraph(g_rho_fin);
	TGraph* g2 = new TGraph(g_j_fin);
	TGraph* g3 = new TGraph(g_psi);

	mg->Add(g1);
	mg->Add(g2);
	mg->Add(g3);
	mg->Draw("apl");
	mg->GetXaxis()->SetTitle("x");
	mg->GetYaxis()->SetTitle("y");
	mg->GetHistogram()->SetTitle("Comparison of probability density with current for a Gaussian (sig=0.1,m=1,t=0)");

	TLegend* legend = new TLegend(0.78, 0.695, 0.98, 0.775);
	legend->AddEntry(g_rho_fin, "Probability density", "l");
	legend->AddEntry(g_j_fin, "Probability current", "l");
	legend->AddEntry(g_psi, "Non-relativistic density", "l");
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
