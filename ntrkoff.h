bool checkBounds(double pt, double eta){
  if( fabs(eta) > 2.4 ){return false;}
  if( pt < 0 || pt > 500 ){return false;}
  return true;
}

double getTrkCorrWeight(TH2 *eff_factor, double pT, double eta){
  if( !checkBounds(pT, eta) ) return 0;
  double factor = 1.0;
  double eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(eta),eff_factor->GetYaxis()->FindBin(pT) );
  factor = (1. / eff); // efficiency and fake rate
  return factor;
}

int get_Ntrkoff(int size, float *eta, float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr){
	int Ntrk_off = 0;
	for(int ii=0; ii<size; ii++){ 
        if(pt[ii] <= 0.4) continue;
		if(fabs(eta[ii]) >= 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] != 1) continue;
		if(fabs(pterr[ii]/pt[ii]) >= 0.1) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 3.0) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= 3.0) continue;
		Ntrk_off=Ntrk_off+1;
	}
	return Ntrk_off;
}

float get_Ntrkcorr(TH2 *eff_histo, int size, float *eta, float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr){
	float Ntrk = 0.0;
	for(int ii=0; ii < size; ii++){ 
        if(pt[ii] <= 0.4) continue;
		if(fabs(eta[ii]) >= 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] != 1) continue;
		if(fabs(pterr[ii]/pt[ii]) >= 0.1) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 3.0) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= 3.0) continue;
		Ntrk = Ntrk + getTrkCorrWeight(eff_histo, pt[ii], eta[ii]);
	}
	return Ntrk;
}

float get_Ntrkcorr_tight(TH2 *eff_histo, int size, float *eta, float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr){
	float Ntrk_tight = 0.0;
	for(int ii=0; ii < size; ii++){ 
        if(pt[ii] <= 0.4) continue;
		if(fabs(eta[ii]) >= 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] != 1) continue;
		if(fabs(pterr[ii]/pt[ii]) >= 0.05) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 2.0) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= 2.0) continue;
		Ntrk_tight = Ntrk_tight + getTrkCorrWeight(eff_histo, pt[ii], eta[ii]);
	}
	return Ntrk_tight;
}

float get_Ntrkcorr_loose(TH2 *eff_histo, int size, float *eta, float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr){
	float Ntrk_loose = 0.0;
	for(int ii=0; ii<size; ii++){ 
        if(pt[ii] <= 0.4) continue;
		if(fabs(eta[ii]) >= 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] != 1) continue;
		if(fabs(pterr[ii]/pt[ii]) >= 0.1) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 5.0) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= 5.0) continue;
		Ntrk_loose = Ntrk_loose + getTrkCorrWeight(eff_histo, pt[ii], eta[ii]);
	}
	return Ntrk_loose;
}