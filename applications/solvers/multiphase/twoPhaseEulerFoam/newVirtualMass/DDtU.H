{
    DDtU1 =
        fvc::ddt(U1)
      + fvc::div(phi, U1)
      - fvc::div(phi)*U1;

    DDtU2 =
        fvc::ddt(U2)
      + fvc::div(phi, U2)
      - fvc::div(phi)*U2;
}
