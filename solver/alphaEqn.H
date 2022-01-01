{
    volScalarField divU(fvc::div(phi));
    
    Pair<tmp<volScalarField>> vDotAlphal =
        mixture->vDotAlphal();
    const volScalarField& vDotcAlphal = vDotAlphal[0]();
    const volScalarField& vDotvAlphal = vDotAlphal[1]();
    const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);

    tmp<surfaceScalarField> talphaPhi;

	fvScalarMatrix alpha1Eqn
	(
		fv::gaussConvectionScheme<scalar>
		(
			mesh,
			phi,
			upwind<scalar>(mesh, phi)
		).fvmDiv(phi, alpha1)
	  - fvm::Sp(divU, alpha1)
	  - fvm::laplacian(alfaCoeffDiffusion_ ,alpha1)
	 ==
		fvm::Sp(vDotvmcAlphal, alpha1)
	  - vDotvAlphal
	);

	alpha1Eqn.solve();

	Info<< "Phase-1 volume fraction = "
		<< alpha1.weightedAverage(mesh.Vsc()).value()
		<< "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
		<< "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
		<< endl;
		
	alpha1Eqn.relax();

	talphaPhi = alpha1Eqn.flux();
    
    alpha2 = 1.0 - alpha1;

    rhoPhi = talphaPhi()*(rho1 - rho2) + phi*rho2;
    
    rho == alpha1*rho1 + alpha2*rho2;

    Info<< "Vapor phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
