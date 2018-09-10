Emp_SPA_pdf_t=function(t,time, status){ 
    #time=observed time vector
    #status=indicate whether observation is censored (0) or uncensored (1)
    #To avoid negative values for log transformation
    Y=log(time+1)
    
    # sort the times, and the status keeping it with the corresp time (concomitant) 
    delta=status[sort.list(Y)]; z=sort(Y); n=length(z);
    #############################################################
    # Compute fhat(Z(i)) when delta(i)=1 :call it as pw[i]      #
    #############################################################
    pw=rep(0,n); pw[1]=delta[1]/n; prod.term=1;
    for (i in 2:n){
        prod.term = prod.term * (1 - delta[i-1]/(n-i+2))
        pw[i] = prod.term * delta[i]/(n-i+1)
    } 
    
    #############################################################
    # Compute phi_n                                             #
    #############################################################
    p.str=sum(pw);
    # decide whether or not to exp tail completetion based on delta[n]
    phi=0; if (delta[n]==0) {phi = -log(1-p.str)/z[n]};
    
    
    #############################################################
    # Empirical MGF of Z                                        #
    #############################################################
    M.mgf=function(s){
        Mds=sum(pw*exp(s*z)); Mcs=0;
        if (delta[n]==0) {Mcs=phi*exp(z[n]*(s-phi))/(phi-s)};
        Ms=Mds+Mcs;
        return(Ms)
    }
    
    
    #############################################################
    # Empirical CGF of Z                                        #
    #############################################################
    K.cgf=function(s){
        Mds=sum(pw*exp(s*z)); Mcs=0;
        if (delta[n]==0) {Mcs=phi*exp(z[n]*(s-phi))/(phi-s)};
        log(Mds+Mcs)
    }
    
    #############################################################
    # 2nd derivative of CGF of Z                                #
    #############################################################
    
    K2.exp=function(s){
        Mds=sum(pw*exp(s*z)); Mcs=0;
        if (delta[n]==0) {Mcs=phi*exp(z[n]*(s-phi))/(phi-s)};
        Ms=Mds+Mcs;
        #   first primes
        Mds.prime=sum(pw*z*exp(s*z)); Mcs.prime=0;
        if (delta[n]==0) {Mcs.prime=Mcs*(z[n]+1/(phi-s))};
        Ms.prime=Mds.prime+Mcs.prime
        #   double primes
        Mds.prime2=sum(pw*z^2*exp(s*z)); Mcs.prime2=0;
        if (delta[n]==0) {Mcs.prime2=Mcs.prime*(z[n]+1/(phi-s))+Mcs/(phi-s)^2};
        Ms.prime2=Mds.prime2+Mcs.prime2
        #   put it all together
        K.prime2=(Ms*Ms.prime2-Ms.prime^2)/Ms^2
        K.prime2
    }
    ##Note: K2.exp==NaN if Ms will be 0. Ms will be zero either if s->-infinity or z-> -infinity
    ##but s should be in the neighborhood of zero. So when optimize CGF adjusted, make a condition that Ms will not be ==0
    
    #############################################################
    # Daniel's saddlepoint pdf at t with s=sh                   #
    #############################################################
    
    fh.com=function(t,s){
        K.exp=K.cgf(s);
        pdf.da=exp(K.exp-t*s)/sqrt(2*pi*K2.exp(s))
        return(pdf.da)
    }
    
    #############################################################
    # Adjusted CGF to optimze (minimize): to find the sh        #
    #############################################################
    
    K.adj.exp=function(s){
        Mds=sum(pw*exp(s*z)); Mcs=0;
        
        if (delta[n]==0){
            Mcs=phi*exp(z[n]*(s-phi))/(phi-s)
        }
        Ms=Mds+Mcs
        
        if(round(Ms,3)==0|round(K2.exp(s),3)<=0){
            return(10000)
        }
        else{
            return(K.cgf(s)-log(t+1)*s)
        }
    }
    
    #############################################################
    # Plot the adjusted CGF                                     #
    #############################################################
    #   s.grid=seq(-5,phi,.01); m=length(s.grid); k.fun=rep(0,m);
    #   for (i in 1:m){k.fun[i]=K.adj.exp(s.grid[i])} 
    #   plot(s.grid,k.fun)
    
    
    #############################################################
    #       Use minimization to find the saddlepoint            #
    #############################################################
    saddle.com<-function(){
        suppressWarnings(optim(0,K.adj.exp,method="Nelder-Mead",hessian=FALSE)$par)
    }
    
    
    sh=saddle.com()
    
    return((fh.com(log(t+1),sh))/(t+1))
    
}
