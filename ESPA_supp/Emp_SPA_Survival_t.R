Emp_SPA_Survival_t=function(t,time, status,inc=.01){ 
#inc=increment given by the user to plot the survival function. 
#This increment is used to avoid the estimate survival at the mean
#To avoid negative values for log transformation

    Y=log(time+1)

#   sort the times, and the status keeping it with the corresp time (concomitant)
    delta=status[sort.list(Y)]; z=sort(Y); n=length(z);
  
#   Compute fhat(Z(i)) when delta(i)=1 :call it as pw[i]
    pw=rep(0,n); pw[1]=delta[1]/n; prod.term=1;
    for (i in 2:n){
        prod.term = prod.term * (1 - delta[i-1]/(n-i+2))
        pw[i] = prod.term * delta[i]/(n-i+1)
    } 
  
#   Compute phi_n    
    p.str=sum(pw);
#   decide whether or not to exp tail completetion based on delta[n]
    phi=0; if (delta[n]==0) {phi = -log(1-p.str)/z[n]};
  
  
#   Empirical MGF of Z 
    M.mgf=function(s){
        Mds=sum(pw*exp(s*z)); Mcs=0;
        if (delta[n]==0) {Mcs=phi*exp(z[n]*(s-phi))/(phi-s)};
        Ms=Mds+Mcs;
        return(Ms)
    }
  

#   Mean of Z (will be used in L & R formula)                 
#   locate mean of emp mgf to get around the little "blip" there         

    M1.prime=function(s){
        Mds=sum(pw*exp(s*z)); Mcs=0;
        if (delta[n]==0) {Mcs=phi*exp(z[n]*(s-phi))/(phi-s)};
        Mds.prime=sum(pw*z*exp(s*z)); Mcs.prime=0;
        if (delta[n]==0) {Mcs.prime=Mcs*(z[n]+1/(phi-s))};  
        M.prime=Mds.prime+Mcs.prime
        return(M.prime)
    }
    emp.mean=M1.prime(0)
  
  
#   Empirical CGF of Z  
    K.cgf=function(s){
        Mds=sum(pw*exp(s*z)); Mcs=0;
        if (delta[n]==0) {Mcs=phi*exp(z[n]*(s-phi))/(phi-s)};
        log(Mds+Mcs)
    }
  
  
#   2nd derivative of CGF of Z 
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
        return(K.prime2)
    }
  
  
#   3rd derivative of CGF of Z  
    K3.exp=function(s){
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
        #   triple primes
        Mds.prime3=sum(pw*z^3*exp(s*z)); Mcs.prime3=0;
        if (delta[n]==0) {
          Mcs.prime3=Mcs.prime2*(z[n]+1/(phi-s))+2*Mcs.prime/(phi-s)^2+2*Mcs/(phi-s)^3
        }
        Ms.prime3=Mds.prime3+Mcs.prime3
        #   put it all together
        K.prime3=(Ms^2*Ms.prime3-3*Ms*Ms.prime*Ms.prime2+2*Ms.prime^3)/Ms^3
        return(K.prime3)
    }
      
 
#   4th derivative of CGF of Z 
    K4.exp=function(s){
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
        # 3rd primes
        Mds.prime3=sum(pw*z^3*exp(s*z)); Mcs.prime3=0;
        if (delta[n]==0) {
          Mcs.prime3=Mcs.prime2*(z[n]+1/(phi-s))+2*Mcs.prime/(phi-s)^2+2*Mcs/(phi-s)^3
        }
        Ms.prime3=Mds.prime3+Mcs.prime3
        # 4th primes
        Mds.prime4=sum(pw*z^4*exp(s*z)); Mcs.prime4=0;
        if (delta[n]==0) {
            Mcs.prime4=Mcs.prime3*(z[n]+1/(phi-s)) + 3*Mcs.prime2/(phi-s)^2
                       + 6*Mcs.prime/(phi-s)^3 + 6*Mcs/(phi-s)^4
        }
        Ms.prime4=Mds.prime4+Mcs.prime4
        #   put it all together
        K.prime4=Ms.prime4/Ms-4*Ms.prime*Ms.prime3/Ms^2+12*Ms.prime2*Ms.prime^2/Ms^3
                 -3*Ms.prime2^2/Ms^2-6*Ms.prime^4/Ms^4
        return(K.prime4)
    }

#   Lugganani & Rice saddlepoint CDF at t with s=sh
    LR.cdf=function(t,s){
        K.exp=K.cgf(s); u=s*sqrt(K2.exp(s)); w=sign(s)*sqrt(2*(s*t-K.exp));
        cdf.here = pnorm(w)+dnorm(w)*(1/w-1/u)
        return(cdf.here)
    }
  
  
#   Use minimization to find the saddlepoint 
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

 
    saddle.me<-function(){
        suppressWarnings(optim(0, K.adj.exp, method="Nelder-Mead", hessian=FALSE)$par)
    }
   
    sh=saddle.me()
    out.list = 1-LR.cdf(log(t+1),sh)
    return(out.list)
}
