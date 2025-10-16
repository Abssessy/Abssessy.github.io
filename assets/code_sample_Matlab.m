function main()
    % Data Description
    % a panel of 48 observations from 1970 to 1986
    % number of observations : 816
    % observation : regional
    % country : United States
    % A data frame containing:
    % column A: observation
    % column B: the state
    % column C: the year
    % column D: gross state products
    % column E: private capital stock
    % column F: public capital
    % column G: labor input measured by the employment in non-agricultural payrolls
    % column H: state unemployment rate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear
    tic
    format short

    disp('Munnel Data: n=48, and T=17 in total')
    data = importdata('Munnel.xls');
    W1   = importdata('weight_Munnel.xls'); % spatial weight matrix
    n    = length(W1);
    it   = 6;     % ending period of the date to be used
    id0  = n*(it-6)+1;  id1 = n*it;
    Y    = data(id0:id1,4);
    X    = data(id0:id1,5:8);

    vnames = char('private capital','public capital','labour input',...
               'unemployment rate   ', 'sigv^2','rho', 'La1', 'La2','La3');

    [nT,k] = size(X);
    n      = length(W1);
    t      = nT/n-1;
    eps    = 0.00001;  % for numerical derivatives of D and D_1.

    nt    = n*(t-1);
    In    = eye(n);
    On    = zeros(n,n);
    It    = eye(t-1);
    Jn    = ones(n,1);
    Jt    = ones(t-1,1);
    Jtt   = ones(t-1,t-1);
    Jnt   = ones(nt,1);
    JtIn  = kron(Jt,In);
    JttIn = kron(Jtt,In);
    gb    = zeros(n,k);

    % Creating the matrix C0 for M-Estimator
    C0  = Fn_CMat(t);
    CC  = kron(C0,In);
    IC0 = pinv(C0);
    IC  = kron(IC0, In);

    W2 = W1;
    W3 = W1;

    WW1 = W1'*W1;
    WT1 = kron(It,W1);
    WT2 = kron(It,W2);
    WT3 = kron(It,W3);

    Y0   = Y(1:n);
    DY   = Y((2*n+1):nT)-Y((n+1):n*t);
    DY_1 = Y((n+1):n*t)-Y(1:(n*(t-1)));
    Dy1  = DY_1(1:n);
    DX   = X((2*n+1):nT,:)-X((n+1):n*t,:);

    %FT   = ones(n,1);
    %%for it = 2:(t-1)
    %    FT = blkdiag(FT,ones(n,1));
    %end
    %DX = [DX FT]; 
    %k  = length(DX(1,:));

    XXC  = [DX, DY_1, WT2*DY_1];

    ev1  = eig(W1);
    ev2  = ev1;
    lb   = 1/min(ev1) + 0.0001;
    ub   = 1/max(ev1) - 0.0001;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the M-Estimators of lambda and rho; then beta and sigma_v^2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    start = [0.01 0.01 0.01 0.01];
    delta = FnAQSE_STLE(DY,DY_1,DX,W1,W2,W3,C0,start); 
    rhoh1 = delta(1);  lamh1 = delta(2); lamh2=delta(3); lamh3=delta(4); 

    B1  = In-lamh1*W1;
    B2  = rhoh1*In+lamh2*W2;
    B3  = In-lamh3*W3;
    IB1 = pinv(B1);
    IB3 = pinv(B3);
    G1  = W1*IB1;
    A3  = W3'*B3+B3'*W3;
    G3  = W3*IB3;
    BB  = IB1*B2;

    DDY   = kron(It,B1)*DY-kron(It,B2)*DY_1;
    IOm   = kron(IC0,B3'*B3);
    beth1 = pinv(DX'*IOm*DX)*DX'*IOm*DDY;
    du    = DDY-DX*beth1;
    sigh1 = sqrt(du'*IOm*du/nt);

    Theta_M  = [beth1',sigh1^2,delta]';

    % Compute Matrices D_1 and D, T >= 2, at rhoh1, lamh1, lamh2
    dpar = [rhoh1, lamh1, lamh2]; 
    D    = Fn_DMat(W1,W2,t,dpar);
    D_1  = Fn_DMat_1(W1,W2,t,dpar);

    % Compute the OPMD Estimate of Var(AQS), T >= 2  
    RM   = BB;
    RM_1 = In;
    for it = 2:(t-1)
        RM   = blkdiag(RM, BB^it);
        RM_1 = blkdiag(RM_1, BB^(it-1));
    end
    BM   = zeros(nt,nt);
    BM_1 = zeros(nt,nt);
    for is = 1:(t-1)
        ids0 = n*(is-1)+1;  ids1 = n*is;  
        BM(ids0:ids1,ids0:ids1) = In;
        for it = (is+1):(t-1)
            idt0 = n*(it-1)+1;  idt1 = n*it;
            BM(idt0:idt1,ids0:ids1) = BB^(it-is);
            BM_1(idt0:idt1,ids0:ids1) = BB^(it-is-1);
        end
    end
    Et   = BM*kron(It,IB1)*DX*beth1;
    Et_1 = BM_1*kron(It,IB1)*DX*beth1;
    t31  = IB1*IB3;
    SM   = BM*kron(It,t31);
    SM_1 = BM_1*kron(It,t31);

    CB   = kron(IC0,B3);
    Pi1  = CB'*DX/sigh1^2;
    Pi2  = CB'*Et_1/sigh1^2;
    Pi3  = CB'*WT1*Et/sigh1^2;
    Pi4  = CB'*WT2*Et_1/sigh1^2;
    Phi1 = kron(IC0,In)/(2*sigh1^4);
    Phi2 = CB*SM_1/sigh1^2;
    Phi3 = CB*WT1*SM/sigh1^2;
    Phi4 = CB*WT2*SM_1/sigh1^2;
    Phi5 = kron(IC0,(G3'+G3))/(2*sigh1^2);
    Psi1 = CB*RM_1/sigh1^2;
    Psi2 = CB*WT1*RM/sigh1^2;
    Psi3 = CB*WT2*RM_1/sigh1^2;
    dv   = kron(It,B3)*du;
    dvv  = reshape(dv,n,t-1);

    dxi1 = zeros(n,t-1);
    dxi2 = zeros(n,t-1);
    dxi3 = zeros(n,t-1);
    dxi4 = zeros(n,t-1);
    dxi5 = zeros(n,t-1);

    for it = 1:(t-1)
       it1  = (it-1)*n+1; it2 = it*n;  
       for is = 1:(t-1)
          is1 = (is-1)*n+1; is2 = is*n;
          P1ts = Phi1(it1:it2, is1:is2);     P1st = Phi1(is1:is2, it1:it2);
          P2ts = Phi2(it1:it2, is1:is2);     P2st = Phi2(is1:is2, it1:it2);
          P3ts = Phi3(it1:it2, is1:is2);     P3st = Phi3(is1:is2, it1:it2);
          P4ts = Phi4(it1:it2, is1:is2);     P4st = Phi4(is1:is2, it1:it2);
          P5ts = Phi5(it1:it2, is1:is2);     P5st = Phi5(is1:is2, it1:it2);
          dvvs = dvv(:,is);
          dxi1(:,it) = dxi1(:,it)+(triu(P1st)'+tril(P1ts)-diag(diag(P1ts)))*dvvs; 
          dxi2(:,it) = dxi2(:,it)+(triu(P2st)'+tril(P2ts)-diag(diag(P2ts)))*dvvs; 
          dxi3(:,it) = dxi3(:,it)+(triu(P3st)'+tril(P3ts)-diag(diag(P3ts)))*dvvs; 
          dxi4(:,it) = dxi4(:,it)+(triu(P4st)'+tril(P4ts)-diag(diag(P4ts)))*dvvs; 
          dxi5(:,it) = dxi5(:,it)+(triu(P5st)'+tril(P5ts)-diag(diag(P5ts)))*dvvs; 
       end
    end

    Pta1 = zeros(n,n);
    Pta2 = zeros(n,n);
    Pta3 = zeros(n,n);
    for is = 1:(t-1)
        is1  = (is-1)*n+1; is2 = is*n;
        Pta1 = Pta1 + Psi1(1:n, is1:is2);
        Pta2 = Pta2 + Psi2(1:n, is1:is2);
        Pta3 = Pta3 + Psi3(1:n, is1:is2);
    end
    Thet1 = Pta1*IB1*IB3;
    Thet2 = Pta2*IB1*IB3;
    Thet3 = Pta3*IB1*IB3;

    gr1c = zeros(n,1);
    g11c = zeros(n,1);
    g21c = zeros(n,1);
    for it = 2:(t-1)
        Pta1 = zeros(n,n);
        Pta2 = zeros(n,n);
        Pta3 = zeros(n,n);
        it1  = (it-1)*n+1; it2 = it*n;  
        for is = 1:(t-1)
            is1  = (is-1)*n+1; is2 = is*n;
            Pta1 = Pta1 + Psi1(it1:it2, is1:is2);
            Pta2 = Pta2 + Psi2(it1:it2, is1:is2);
            Pta3 = Pta3 + Psi3(it1:it2, is1:is2);
        end 
        gr1c = gr1c + dvv(:,it).*(Pta1*Dy1);
        g11c = g11c + dvv(:,it).*(Pta2*Dy1);
        g21c = g21c + dvv(:,it).*(Pta3*Dy1);
    end

    for j=1:k
        gb(:,j) = sum(reshape(Pi1(:,j),n,t-1).*dvv, 2); 
    end
    gs  = sum(dvv.*dxi1, 2)-(t-1)/(2*sigh1^2); 
    gl3 = sum(dvv.*dxi5, 2)-(t-1)*diag(G3);

    Dy1o = B3*B1*Dy1;
    dvv1 = dvv(:,1);
    dr   = sigh1^2*sum(reshape(diag(CC*Phi2),n,t-1), 2);
    dl1  = sigh1^2*sum(reshape(diag(CC*Phi3),n,t-1), 2);
    dl2  = sigh1^2*sum(reshape(diag(CC*Phi4),n,t-1), 2);

    tem  = triu(Thet1)'+tril(Thet1)-2*diag(diag(Thet1));
    gr1a = dvv1.*(tem*Dy1o);
    gr1b = diag(Thet1).*(dvv1.*Dy1o+sigh1^2);
    gr2  = sum(reshape(Pi2,n,t-1).*dvv, 2);
    gr3  = sum(dvv.*dxi2, 2) - dr; 
    gr   = gr1a + gr1b + gr1c + gr2 + gr3;

    tem = triu(Thet2)'+tril(Thet2)-2*diag(diag(Thet2));
    g11a= dvv1.*(tem*Dy1o);
    g11b= diag(Thet2).*(dvv1.*Dy1o+sigh1^2);
    g12 = sum(reshape(Pi3,n,t-1).*dvv, 2);
    g13 = sum(dvv.*dxi3, 2) - dl1; 
    gl1 = g11a + g11b + g11c + g12 + g13;

    tem = triu(Thet3)'+tril(Thet3)-2*diag(diag(Thet3));
    g21a= dvv1.*(tem*Dy1o);
    g21b= diag(Thet3).*(dvv1.*Dy1o+sigh1^2);
    g22 = sum(reshape(Pi4,n,t-1).*dvv, 2);
    g23 = sum(dvv.*dxi4, 2) - dl2; 
    gl2 = g21a + g21b + g21c + g22 + g23;

    gv  = [gb, gs, gr, gl1, gl2, gl3];

    OPMD = gv'*gv;
    VSC  = pinv(OPMD);

    % Compute the Hessian Matrix of the AQS Function    
    % Derivative of D_1 w.r.t. rho:
    if t == 2
        D_1r = On;
    elseif t == 3
        D_1r = [On, On; IB1^2, On];
    else  %using numerical derivatives
        dpar = [rhoh1+eps,lamh1,lamh2];
        D_1r = (Fn_DMat_1(W1,W2,t,dpar) - D_1)/eps;
    end

    % Derivative of D_1 and D w.r.t. lambda1, lambda2:
    IWI1 = IB1*W1*IB1;
    IWI2 = IB1*W2*IB1;
    te1  = IWI1*B2*IB1+BB*IWI1-2*IWI1;
    te2  = (IB1*W1*BB^2+BB*IB1*W1*BB-2*IB1*W1*BB)*IB1;
    te3  = (IB1*W2*BB+BB*IB1*W2-2*IB1*W2)*IB1; 
    if t == 2
        D_11 = IWI1;    D_12 = On;
        DL1  = te1;     DL2  = IB1*W2*IB1;
    elseif t == 3    
        D_11 = [IWI1, On; te1, IWI1];
        D_12 = [On, On; IWI2, On];
        %DL1  = [te1, IWI1; (IB1*W1*BB^2+BB*IB1*W1*BB-2*IB1*W1*BB)*IB1, te1];
        DL1  = [te1, IWI1; te2+(In-BB)^2*IWI1, te1];
        DL2  = [IWI2, On; te3, IWI2];
    else  % using numerical derivatives
        dpar = [rhoh1,lamh1+eps,lamh2];
        DL1  = (Fn_DMat(W1,W2,t,dpar)-D)/eps;
        D_11 = (Fn_DMat_1(W1,W2,t,dpar)-D_1)/eps;

        dpar = [rhoh1,lamh1,lamh2+eps];
        DL2  = (Fn_DMat(W1,W2,t,dpar)-D)/eps;
        D_12 = (Fn_DMat_1(W1,W2,t,dpar)-D_1)/eps;
    end

    dIOm = kron(IC0,A3);
    Hbb  = -DX'*IOm*DX/sigh1^2;
    Hbs  = -DX'*IOm*du/sigh1^4;
    Hbr  = -DX'*IOm*DY_1/sigh1^2;
    Hbl1 = -DX'*IOm*WT1*DY/sigh1^2;
    Hbl2 = -DX'*IOm*WT2*DY_1/sigh1^2;
    Hbl3 = -DX'*dIOm*du/sigh1^2;

    Hss  = -du'*IOm*du/sigh1^6 + 0.5*nt/sigh1^4; 
    Hsr  = -DY_1'*IOm*du/sigh1^4;
    Hsl1 = -DY'*WT1'*IOm*du/sigh1^4;
    Hsl2 = -DY_1'*WT2'*IOm*du/sigh1^4;
    Hsl3 = -0.5*du'*dIOm*du/sigh1^4;

    Hrr  = -DY_1'*IOm*DY_1/sigh1^2 + trace(IC*D_1r);
    Hrl1 = -DY_1'*IOm*WT1*DY/sigh1^2 + trace(IC*D_11);
    Hrl2 = -DY_1'*IOm*WT2*DY_1/sigh1^2 + trace(IC*D_12);
    Hrl3 = -DY_1'*dIOm*du/sigh1^2;

    Hl11 = -DY'*WT1'*IOm*WT1*DY/sigh1^2 + trace(IC*DL1*WT1);          
    Hl12 = -DY'*WT1'*IOm*WT2*DY_1/sigh1^2 + trace(IC*DL2*WT1);          
    Hl13 = -DY'*WT1'*dIOm*du/sigh1^2;  

    Hl22 = -DY_1'*WT2'*IOm*WT2*DY_1/sigh1^2 + trace(IC*D_12*WT2);
    Hl23 = -DY_1'*WT2'*dIOm*du/sigh1^2;
    Hl33 = -du'*kron(IC0, W3'*W3)*du/sigh1^2 - (t-1)*trace(G3^2);

    HH   = [Hbb,   Hbs,  Hbr,  Hbl1, Hbl2, Hbl3; 
            Hbs',  Hss,  Hsr,  Hsl1, Hsl2, Hsl3; 
            Hbr',  Hsr,  Hrr,  Hrl1, Hrl2, Hrl3; 
            Hbl1', Hsl1, Hrl1, Hl11, Hl12, Hl13;
            Hbl2', Hsl2, Hrl2, Hl12, Hl22, Hl23;
            Hbl3', Hsl3, Hrl3, Hl13, Hl23, Hl33]; 

    SIGH = -HH;   IHH = pinv(SIGH);    dIHH = diag(IHH);

    se_M  = sqrt(diag(IHH*OPMD*IHH));
    t_M  = Theta_M./se_M;

    result = [Theta_M se_M t_M];
    disp('Estimation and inference for Munnel Panel Data, STLE model')
    disp('Response: Gross social product')
    disp(' ')
    disp('                     M-Est     se_M      t_M')
    disp([vnames, num2str(result, '% 10.4f')])
    toc
end

function f = Fn_CMat(t)
    % For evaluating the matrices C0 and C(varsigma), here t represents T (>=3)
    C0 = zeros(t-1,t-1);
    C0(1,1) = 2;  C0(1,2) = -1;
    C0(t-1,t-2) = -1; C0(t-1,t-1) = 2;
    for it = 2:(t-2)
        C0(it,it) = 2;
        C0(it,it-1) = -1; 
        C0(it,it+1) = -1;
    end
    f = C0;
end

function f = Fn_DMat(W1,W2,t,dpar)
    % Compute Matrices D, T >= 2, at rhoh1, lamh1, lamh2
    n   = length(W1);
    nt  = n*(t-1);
    In  = eye(n);
    B1  = In-dpar(2)*W1;
    B2  = dpar(1)*In + dpar(3)*W2;
    IB1 = pinv(B1);
    BB  = IB1*B2;
    te0 = (BB-2*In)*IB1;
    te1 = (In-BB)^2*IB1;

    D  = zeros(nt,nt);
    if t == 2
        D  = te0;
    elseif t == 3
        D  = [te0, IB1; te1, te0];
    else
        for is = 1:(t-2);
            ids0 = n*(is-1)+1;  ids1 = n*is;  
            D(ids0:ids1,ids0:ids1) = te0;
            D(ids0:ids1,(ids1+1):(ids1+n)) = IB1;
            D((ids1+1):(ids1+n),ids0:ids1) = te1;
            for it = (is+2):(t-1);
                idt0 = n*(it-1)+1;  idt1 = n*it;
                D(idt0:idt1,ids0:ids1) = BB^(it-is-1)*te1;
            end
        end
        D((n*(t-2)+1):n*(t-1), (n*(t-2)+1):n*(t-1)) = te0;
    end

    f = D;
end

function f = Fn_DMat_1(W1,W2,t,dpar)
    % Compute Matrices D, T >= 2, at rhoh1, lamh1, lamh2
    n   = length(W1);
    nt  = n*(t-1);
    In  = eye(n);
    On  = zeros(n,n);
    B1  = In-dpar(2)*W1;
    B2  = dpar(1)*In + dpar(3)*W2;
    IB1 = pinv(B1);
    BB  = IB1*B2;
    te0 = (BB-2*In)*IB1;
    te1 = (In-BB)^2*IB1;

    D_1 = zeros(nt,nt);
    if t == 2
        D_1 = IB1;
    elseif t == 3
        D_1 = [IB1, On; te0, IB1];
    else
        for is = 1:(t-2);
            ids0 = n*(is-1)+1;  ids1 = n*is;  
            D_1(ids0:ids1,ids0:ids1) = IB1;
            D_1((ids1+1):(ids1+n),ids0:ids1) = te0;
            for it = (is+2):(t-1);
                idt0 = n*(it-1)+1;  idt1 = n*it;
                D_1(idt0:idt1,ids0:ids1) = BB^(it-is-2)*te1;
            end
        end
        D_1((n*(t-2)+1):n*(t-1), (n*(t-2)+1):n*(t-1)) = IB1;
    end

    f = D_1;
end

function [x fval] = FnAQSE_STLE(DY,DY_1,DX,W1,W2,W3,C0,par)
    % Nested function that computes the objective function
    warning('off', 'all');
    opd = optimset('Display', 'off');

    [x fval] = fsolve(@ScSDPD_F, par, opd);

    function f = ScSDPD_F(par)
        n  = length(W1);
        t  = 1+length(DY)/n;
        In = eye(n);
        It = eye(t-1);
        rho = par(1);
        la1 = par(2);
        la2 = par(3);
        la3 = par(4);
        WT1 = kron(It,W1);
        WT2 = kron(It,W2);
        B1  = In-la1*W1;
        B2  = rho*In+la2*W2;
        B3  = In-la3*W3;
        IC0 = pinv(C0);
        IC  = kron(IC0,In);
        IB3 = pinv(B3);
        IOm = kron(IC0,B3'*B3);
        DDY = kron(It,B1)*DY-kron(It,B2)*DY_1;
        du  = DDY-DX*pinv(DX'*IOm*DX)*DX'*IOm*DDY;
        sigh = sqrt(du'*IOm*du/(n*(t-1)));

        % Compute Matrices D_1 and D, T >= 2, at rhoh1, lamh1, lamh2
        dpar = [rho, la1, la2];
        D    = Fn_DMat(W1,W2,t,dpar);
        D_1  = Fn_DMat_1(W1,W2,t,dpar);

        f(1) = du'*IOm*DY_1/sigh^2 + trace(D_1*IC);
        f(2) = du'*IOm*WT1*DY/sigh^2 + trace(D*IC*WT1);
        f(3) = du'*IOm*WT2*DY_1/sigh^2 + trace(D_1*IC*WT2);
        f(4) = du'*kron(IC0, W3'*B3)*du/sigh^2 - (t-1)*trace(W3*IB3);
    end
end