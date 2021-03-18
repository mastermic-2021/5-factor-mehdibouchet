MAX_BOUND= 2^15;
INIT_BOUND= 2^5;

encode(m)=fromdigits(Vec(Vecsmall(m)),128);
decode(m)=Strchr(digits(m, 128));
[n,c] = readvec("input.txt");


\\ **********************************************************************************
\\  Fonctions de l'exercice :                                                       *
\\                                                                                  *
\\    dechiffre(c, N): Déchiffrement d'un <chiffré> c par le système de Rabin       *
\\                                                                                  *
\\    rho(n): Méthode p-1 de Pollard. Si n=p*q et p-1 B-superfriable                *
\\      en posant M = B!, p-1 | M donc                                              *
\\      x= a^M mod p = a ^ (p-1) mod p = 1 mod p.                                   *
\\      Donc  p | a^M - 1 et p | p*q=n donc pgcd(a^M - 1, n) = p                    *
\\                                                                                  *
\\  Méthode:                                                                        *
\\    On va bruteforce la borne de friabilité B.                                    *
\\                                                                                  *
\\ **********************************************************************************


rho(n) = {
    my(B, M, p, a);
    B= INIT_BOUND;
    a= Mod(3, n);

    while(B < MAX_BOUND,
        M= B!;
        x= a^M;
        p= gcd(lift(x)-1,n);        
        if( p > 1, 
            if( p == n, a= Mod(random(n), n); B= INIT_BOUND, return(p)) );
        B*= 2;
    );
    -1;
}

dechiffre(c, N)= {
    my(n);
    [c, b1, b2]= c;
    [p, q]= N;
    n= p*q;

    c= Mod(c, n);
    m_p= c ^( (p + 1) / 4 );
    m_q= c ^( (q + 1) / 4 );
    [u,v]= bezout(p,q);

    if( b2 == 1,
        r= u*p*m_q + v*q*m_p;
        if( b1 == 1, return( lift(r) ), return( lift(n-r) ) ),
        s= u*p*m_q - v*q*m_p;
        if( b1 == 1, return( lift(s) ), return( lift(n-s) ) );
    );
}

p= rho(n);
if( p == -1, print("Rho n'a pas fonctionné..."); quit(-1) );

q= n/p;
m= decode( dechiffre(c, [p,q]) );
print( m );