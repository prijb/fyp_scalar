#include <vector>

double int_rho(double *xvals, double *params){
	
	double pos = xvals[0];
	double t = params[0];
        double sig_p = params[1];
        double m = params[2];
	double pi = TMath::Pi();

	TF2* rho = new TF2("rho", "(0.25/pi)*(sqrt(2./pi)*(1/[1]))*(1 + [3]*[3]/(sqrt([3]*[3] + x*x)*sqrt([3]*[3] + y*y)) + x*y/(sqrt([3]*[3] + x*x)*sqrt([3]*[3] + y*y)))*cos((x-y)*[0] - (sqrt([3]*[3] + x*x) - sqrt([3]*[3] + y*y))*[2])*TMath::Exp(-1.*(x*x + y*y)/([1]*[1])) ", -100*sig_p, 100*sig_p, -100*sig_p, 100*sig_p);
	rho->SetParameter(0,pos);
	rho->SetParameter(1,sig_p);
	rho->SetParameter(2,t);
	rho->SetParameter(3,m);
	
	int  n_points = params[3];
	double px_start = -10*sig_p;
	double px_end = 10*sig_p;
	double py_start = -10*sig_p;
        double py_end = 10*sig_p;
	double dp = (px_end-px_start)/n_points;
	double px_a = px_start;
	double px_b = 0;
	double py_a = py_start;
        double py_b = 0;
	double f_a = rho->Eval(px_a,py_a);
	double f_mid = 0;
	double f_b = 0;
	double d_intg = 0;
	double intg = 0;

	for (int i=0; i< n_points; i++){
		for (int k=0; k< n_points; k++){
			py_b = py_a + dp; 
		     	f_b = rho->Eval(px_a,py_b);
			f_mid = rho->Eval(px_a,0.5*(py_a+py_b));
			
			d_intg =(1./6.)*(f_a + 4*f_mid + f_b);

			intg += d_intg;
			
			py_a = py_b;
			f_a = f_b;
		}
		
		px_a = px_a + dp;
		py_a = py_start;
		f_a = rho->Eval(px_b,py_a);
		/* Debugging
                cout << "Iteration: " << i << " p_a: "<< p_a << " p_b: "<< p_b << endl;
                cout << " f_a: "<< f_a << " f_b: "<< f_b << endl;
		cout << "Integral addition: " << d_intg << endl;
                */
	}	

	return dp*dp*intg;
}

double int_j(double *xvals, double *params){

	double pos = xvals[0];
        double t = params[0];
        double sig_p = params[1];
        double m = params[2];
        double pi = TMath::Pi();

        TF2* j = new TF2("j", "(0.25/pi)*(sqrt(2./pi)*(1/[1]))*(x/sqrt([3]*[3] + x*x) + y/sqrt([3]*[3] + y*y))*cos((x-y)*[0] - (sqrt([3]*[3] + x*x) - sqrt([3]*[3] + y*y))*[2])*TMath::Exp(-1.*(x*x + y*y)/([1]*[1])) ", -100*sig_p, 100*sig_p, -100*sig_p, 100*sig_p);
        j->SetParameter(0,pos);
        j->SetParameter(1,sig_p);
        j->SetParameter(2,t);
        j->SetParameter(3,m);

        int  n_points = params[3];
        double px_start = -10*sig_p;
        double px_end = 10*sig_p;
        double py_start = -10*sig_p;
        double py_end = 10*sig_p;
        double dp = (px_end-px_start)/n_points;
        double px_a = px_start;
        double px_b = 0;
        double py_a = py_start;
        double py_b = 0;
        double f_a = j->Eval(px_a,py_a);
        double f_mid = 0;
        double f_b = 0;
        double d_intg = 0;
        double intg = 0;

        for (int i=0; i< n_points; i++){
                for (int k=0; k< n_points; k++){
                        py_b = py_a + dp;
                        f_b = j->Eval(px_a,py_b);
                        f_mid = j->Eval(px_a,0.5*(py_a+py_b));

                        d_intg =(1./6.)*(f_a + 4*f_mid + f_b);

                        intg += d_intg;

                        py_a = py_b;
                        f_a = f_b;
                }

                px_a = px_a + dp;
                py_a = py_start;
                f_a = j->Eval(px_b,py_a);
                /* Debugging
                cout << "Iteration: " << i << " p_a: "<< p_a << " p_b: "<< p_b << endl;
                cout << " f_a: "<< f_a << " f_b: "<< f_b << endl;
                cout << "Integral addition: " << d_intg << endl;
                */
        }

        return dp*dp*intg;
}

void graph_int(){

	double t = 10.;
	double sig_p = 1.;
	double sig_x = 1.;
	//double sig_x = 1./sig_p;
	double m = 1;
	int n_points = 1000;

	double lim = 20*sig_x;
	
	TF1* p_dens = new TF1("p_dens",int_rho,-lim,lim,4);
        TF1* j_dens = new TF1("j_dens",int_j,-lim,lim,4);

	p_dens->SetParameters(t,sig_p,m,n_points);
	j_dens->SetParameters(t,sig_p,m,n_points);

	j_dens->SetLineColor(1);

	double p_int = p_dens->Integral(-lim,lim,1.e-3);
	double j_int = j_dens->Integral(-lim,lim,1.e-3);
	
	std::ostringstream title;
	title << std::fixed << std::setprecision(2) << "Comparison of probability density with current for a Gaussian (sig=" << 1/sig_p << ",m=" << m << ",t=" << t << ")";
	std::string title_graph = title.str();	
	
	cout << "The integral of probability density is:  " << p_int << endl;

	
	TCanvas* c1 = new TCanvas();
        c1->SetGrid();

        TMultiGraph* mg = new TMultiGraph();
        TGraph* g1 = new TGraph(p_dens);
        TGraph* g2 = new TGraph(j_dens);

        mg->Add(g1);
        mg->Add(g2);
        mg->Draw("apl");
        mg->GetXaxis()->SetTitle("x");
        mg->GetYaxis()->SetTitle("y");
        //mg->GetHistogram()->SetTitle("Comparison of probability density with current for a Gaussian (sig=0.1,m=1,t=5)");
	mg->GetHistogram()->SetTitle(title_graph.c_str());

        TLegend* legend = new TLegend(0.78, 0.695, 0.98, 0.775);
        legend->AddEntry(p_dens, "Probability density", "l");
        legend->AddEntry(j_dens, "Probability current", "l");
        legend->Draw("Same");

}	
