
namespace stp_funct {


    // to calculate the effective charge of a particle in
    // stopping power calculations (Spar-Armstrongs &
    // Chandler-ORNL-4869 [1973] )
    double zeff( double z, double beta) {

        double zp = TMath::Power( z, 0.666667 ); // z^(2/3)

        double zeff;

        if( beta <= 0.07*zp ) {

            // low beta case

            zeff = z * ( 1 - TMath::Exp( -125.0*beta/zp ) );

        } else {

            // high beta case.

            zeff = z;
        }


        return zeff;
    }


    double btoep( double dbsq ){

        // output: the energy of a proton [MeV/c^2]
        // input:
        // dbsq ==> beta^2

        double d = TMath::Sqrt( 1/(1-dbsq) ) -1 ;

        double btoep = 938.2592 * d ;

        return btoep;
    }



    double beta2( double E_k, double m ) {

        // output: beta^2
        // input :
        // m = mass  [MeV/c^2]
        // E_total = E_k + mc2.

        double r = E_k/m + 1;
        double beta2 = 1-( 1/r/r);

        return beta2;
    }


    double algip( double z ) {


        double pot[ 13 + 1] = { 0,
                                18.7, 42.0, 39.0,  60.0,  68.0,
                                78.0, 99.5, 98.5, 117.0, 140.0,
                                150.0, 157.0, 163.0 };

        int iz = int( z + 0.05 );


        double potl;

        if( iz > 13 )
        {
            potl = 9.76*z + 58.5/ TMath::Power( z, 0.19 );
        }
        else
        {
            potl = pot[ iz ];
        }

        return TMath::Log( potl * 1.0E-6 );

    }

}


class stp_class {

private:

    double target_density; // [g/cm3]
    double target_A_avg;
    double target_Z_avg;

    double beam_A;
    double beam_Z;
    double beam_energy; // [MeV]

    double dedxmg;  // stoping power, [ MeV/ (mg/cm^2) ]
    double rangemg; // stopping range [ (mg/cm^2) ]

    double avip;
    double avz;
    double elni;
    double eden;


    // for calculation ============
    double dedx( double emass, double epart, double zpart );
    double dedxp( double energy, double db2, double beta );
    double deff( double e);
    double shell( double e);
    double range( double dx,double emass,double epart,double zpart,double sp );
    double de( double dx, double emass, double epart, double zpart,double sp );
    //==============================

public:

    // class constructor
    stp_class(  double in_tgt_d,
                double in_tgt_A, double in_tgt_Z,
                double in_bm_A,  double in_bm_Z,
                double in_bm_eng )
    {
        target_density = in_tgt_d;
        target_A_avg   = in_tgt_A;
        target_Z_avg   = in_tgt_Z;

        beam_A = in_bm_A;
        beam_Z = in_bm_Z;
        beam_energy = in_bm_eng;

        process();
    }


    void   process();  // to do the calculation for dE/dx and range.


    void set_target_density( double density_in_gcm3 ) { target_density = density_in_gcm3; }
    void set_target_A_avg( double A ) { target_A_avg = A; }
    void set_target_Z_avg( double Z ) { target_Z_avg = Z; }
    void set_beam_energy( double eng_in_MeV ) { beam_energy = eng_in_MeV; }
    void set_beam_A( double A ) { beam_A = A; }
    void set_beam_Z( double Z ) { beam_Z = Z; }

    double get_dEdx() { return dedxmg; }
    double get_range() { return rangemg; }

};

double stp_class::de( double dx, double emass, double epart, double zpart,double sp ) {

    // to calculate energy lost over dx cms assuming a quadratic
    //  relationship between e & x.

    double deltae = sp * dx;

    double enew = epart - deltae;

    double de;

    if( enew < 0.0 )
    {
        de = epart;
    }
    else
    {
        double spges = dedx( emass, enew, zpart );

        de = dx * (0.75 * sp + ( 0.25 * spges/ sp ) *spges );
    }

    return de;

}

double stp_class::range( double dx,double emass,double epart,double zpart,double sp ){

    double spe = sp;

    double ed = epart;

    double r = 0.0;

    double es = epart;

    double dxt = dx;

    double dxs = dxt;

    double f = 1.0;

    double g = 1.0;

    double er1;

    double er2 = 1.0;

    int jf;

    double e;

    ed = ed - de( dxt, emass, es, zpart, spe );

    e = ed;

    double range;

    while( 1 ) {

        if( e <=0.0000001 ) { break; }

        r = r + dxt;

        es = e;

        spe = dedx( emass, es, zpart );

        er1 = epart/e;



        if( TMath::Abs( er1/er2 -1.0) < 0.01505 ) {

            g = g *1.1;

            dxt = g * dxs;

            dxs = dxt;

        } else {

            jf = TMath::Min( int(er1), 32  );

            f = 1.0/jf ;
        }

        er2 = er1;

        dxt = f * dxs ;

        ed = ed - de( dxt, emass, es, zpart, spe );

        e = ed;

    } //------------------end of while

    range = r + es / spe;

    return range;
}

double stp_class::shell( double e) {

    // to calculate the shell correction term in the stopping power.

    const double emp = 938.2592;
    const double emp1 = 0.00106580356;
    const double a1 = 0.422377E-6;
    const double b1 = 0.304043E-7;
    const double c1 = -0.38106E-9;

    const double a2 = 3.858019E-9;
    const double b2 = -0.1667989E-9;
    const double c2 = 0.157955E-11;


    const double p1 = 4.774248E-4;
    const double p2 = 1.143478E-3;
    const double p3 = -5.63392E-2;
    const double p4 = 4.763953E-1;
    const double p5 = 4.844536E-1;

    const double w1 = -1.819954E-6;
    const double w2 = -2.232760E-5;
    const double w3 = 1.219912E-4;
    const double w4 = 1.837873E-3;
    const double w5 = -4.457574E-3;
    const double w6 = -6.837103E-2;
    const double w7 = 5.266586E-1;
    const double w8 = 3.743715E-1;


    const double const1 = 0.021769515;
    const double zal = 12.95;
    const double zh2o = 3.34;
    double shell;


    if( e >= 2008.0 ) {

        double gnu2 = 1/( (e*emp1) *(e*emp1+2)  );

        double f1 = gnu2 * ( a1 + gnu2 * (b1 + c1 * gnu2 ) );

        double f2 = gnu2 * ( a2 + gnu2 * (b2 + c2 * gnu2 ) );

        shell = avip * avip * ( f1 + f2* avip )/ avz;

        return shell;
    }

    double be2 = stp_funct::beta2( e, emp );

    shell = const1 + TMath::Log( be2) - TMath::Log( avip );

    double x = 18769.0 * be2 / avz;

    double xlog = TMath::Log( x );

    double xl;

    double xl1, xl2;


    if( avz > zal ) {

        xl = p5 + xlog * ( p4 + xlog * ( p3 + xlog * ( p2 + xlog * p1 ) ) );

        xl = TMath::Exp( xl );

        shell = shell - xl;

    } else {


        xl1 = w8 + xlog * (w7+xlog  *( w6+xlog * ( w5 + xlog * (w4 + xlog * ( w3+ xlog * ( w2 + xlog * w1)))))) ;

        xl1 = TMath::Exp( xl1 );

        if( avz > zh2o ) {

            xl = p5 + xlog * ( p4 + xlog * ( p3 + xlog * ( p2 + xlog * p1 ) ) );

            xl = TMath::Exp( xl );

            xl2 = xl;

            xl = xl1 + (avz-zh2o )/ (zal-zh2o) * ( xl2-xl1);
        }
        else {

            xl = xl1;

        }

        shell = shell - xl;
    }

    return shell;


}

double stp_class::deff( double e) {

    // calculate the density effect correction term in stopping power

    double del = e * 0.00106580356;

    del = TMath::Log( del ) + TMath::Log( del + 2 );

    del = log( 1.378E-9 * eden ) + del - 2.0 * elni - 1.0;


    if( del < 0 ) { return 0; }
    else { return del; }

}


double stp_class::dedxp( double energy, double db2, double beta ) {

    // only valid for beta > 0.0046: 10 KeV protons.

    double dsp = eden * 0.5099147;

    double ze = stp_funct::zeff( 1.0, beta );

    dsp = dsp * (ze*ze)/ db2;

    double d = TMath::Log( 1.022008*db2/(1.0-db2) ) - db2 ;

    double delta = deff( energy );

    double coz = shell ( energy );

    dsp = dsp * ( d - elni - coz - delta*0.5 );

    return dsp;

}

double stp_class::dedx( double emass, double epart, double zpart  ) {

    // only vaild for beta > 0.0046 * z ** (1/3)

    double db2 = stp_funct::beta2( epart, emass);

    double pe = stp_funct::btoep( db2 );

    double beta = TMath::Sqrt( db2 );

    double gog = stp_funct::zeff( zpart, beta) / stp_funct::zeff(1, beta );

    gog = gog * gog;

    double dedx = gog * dedxp( pe, db2, beta );

    return dedx ;

}

void stp_class::process() {



    double num = target_density / ( target_A_avg * 1.660543 );

    double num_Z = num * target_Z_avg;

    double num_Zi = num_Z * stp_funct::algip( target_Z_avg );


    avz = num_Z/num;

    eden = num_Z;

    elni = num_Zi/ num_Z;
    avip = TMath::Exp( elni );

    double tgtionpot = avip * 1000000;

    double emass = beam_A * 931.4812; // in [ MeV/c^2]

    double dx = 0.05/ dedx( emass, beam_energy, beam_Z);

    double sp = dedx( emass, beam_energy, beam_Z);



    double r = range( dx, emass, beam_energy, beam_Z, sp );

    dedxmg = sp*0.001/ target_density;

    rangemg = r * target_density * 1000;


}



/*
  the demo of how to use this class.
*/
void stp() {

    stp_class* demo_obj;

    // target = CH2
    // C ==> A = 12, Z = 6
    // H ==> A =  1, Z = 1
    // A_avg = ( 12*1 + 1*2) / 3
    // Z_avg = (  6*1 + 1*2) / 3
    // beam = 86Kr, A = 86, Z = 36
    
    double target_density = 0.93; // [g/cm3]
    double target_A_avg = 14./3;
    double target_Z_avg = 8./3;
    double beam_A = 86;
    double beam_Z = 36;
    double beam_eng = 90; // [MeV]
    
    demo_obj = new stp_class( target_density, 
                              target_A_avg, 
                              target_Z_avg,
                              beam_A, 
                              beam_Z, 
                              beam_eng );


    // ====== demo1 ======

    cout << "dE/dx = " << demo_obj->get_dEdx() << endl;


    // ====== demo2 ======
    // to show dE/dx as a function of 86Kr energy.

    // use TGraph for graph
    int pointN = 0; 
    TGraph* gr = new TGraph();

    // we would like to wrie out dE/dx data.
    fstream fout;
    fout.open( "86Kr_in_CH2.dat", ios::out );

    for( double eng = 1; eng < 5000; eng+=1) 
    {

        // get dE/dx 
        demo_obj->set_beam_energy( eng );
        demo_obj->process();
        double dEdx = demo_obj->get_dEdx();

        fout << eng << "\t" << dEdx << endl;

        gr->SetPoint( pointN, eng, dEdx );
        pointN+=1;
    }

    gr->Draw("APL");
    fout.close();
    
}
