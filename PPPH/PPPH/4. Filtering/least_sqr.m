function [ xk,pk,kof ] = least_sqr(sno,sls,pno,meas,sit,bp,freq,options)


mno = 2*sno;

x1k = zeros(pno,1);

Hk = zeros(mno,pno);
Zk = zeros(mno,  1);
Ck = zeros(mno,  1);

Rk = eye(mno);

while 1
    for k=1:sno
        sat = meas((4*k)-3,8:10);
        
        rho = norm(sat - x1k(1:3,1)');
        
        Jx = (x1k(1,1) - sat(1,1))/rho;
        Jy = (x1k(2,1) - sat(1,2))/rho;
        Jz = (x1k(3,1) - sat(1,3))/rho;
        Jt = 1;
        Jw = meas((4*k - 3),28);
        Jtn= meas((4*k - 3),29);
        Jte= meas((4*k - 3),30);
        Jn = 1;
        Jr = 1;
        Je = 1;
        Jc = 1;
        
        s = (2*k - 1);
        f = (2*k);
        Hk(s:f,1     ) = Jx;
        Hk(s:f,2     ) = Jy;
        Hk(s:f,3     ) = Jz;
        Hk(s:f,4     ) = Jt;
        Hk(s:f,5     ) = Jw;
        Hk(f  ,(bp+k)) = Jn;
        switch options.TroGrad
            case 0
                if sit==3 || sit==4 || sit==5
                    if sls(k)<33
                        Hk(s:f,6) = 0;
                    else
                        Hk(s:f,6) = 1;
                    end
                elseif sit==6 || sit==7
                    if sls(k)<59
                        Hk(s:f,6) = 0;
                    else
                        Hk(s:f,6) = 1;
                    end
                elseif sit==8 || sit==9
                    if sls(k)<33
                        Hk(s:f,6) = 0;
                        Hk(s:f,7) = 0;
                    elseif sls(k)<59
                        Hk(s:f,6) = 1;
                        Hk(s:f,7) = 0;
                    else
                        Hk(s:f,6) = 0;
                        Hk(s:f,7) = 1;
                    end
                elseif sit==10
                    if sls(k)<33
                        Hk(s:f,6) = 0;
                        Hk(s:f,7) = 0;
                    elseif sls(k)<89
                        Hk(s:f,6) = 1;
                        Hk(s:f,7) = 0;
                    else
                        Hk(s:f,6) = 0;
                        Hk(s:f,7) = 1;
                    end
                elseif sit==11
                    if sls(k)<59
                        Hk(s:f,6) = 0;
                        Hk(s:f,7) = 0;
                    elseif sls(k)<89
                        Hk(s:f,6) = 1;
                        Hk(s:f,7) = 0;
                    else
                        Hk(s:f,6) = 0;
                        Hk(s:f,7) = 1;
                    end
                elseif sit==12
                    if sls(k)<33
                        Hk(s:f,6) = 0;
                        Hk(s:f,7) = 0;
                        Hk(s:f,8) = 0;
                    elseif sls(k)<59
                        Hk(s:f,6) = 1;
                        Hk(s:f,7) = 0;
                        Hk(s:f,8) = 0;
                    elseif sls(k)<89
                        Hk(s:f,6) = 0;
                        Hk(s:f,7) = 1;
                        Hk(s:f,8) = 0;
                    elseif sls(k)<106
                        Hk(s:f,6) = 0;
                        Hk(s:f,7) = 0;
                        Hk(s:f,8) = 1;
                    end
                end
            case 1
                Hk(s:f,6     ) = Jtn;
                Hk(s:f,7     ) = Jte;
                if sit==3 || sit==4 || sit==5
                    if sls(k)<33
                        Hk(s:f,8) = 0;
                    else
                        Hk(s:f,8) = 1;
                    end
                elseif sit==6 || sit==7
                    if sls(k)<59
                        Hk(s:f,8) = 0;
                    else
                        Hk(s:f,8) = 1;
                    end
                elseif sit==8 || sit==9
                    if sls(k)<33
                        Hk(s:f,8) = 0;
                        Hk(s:f,9) = 0;
                    elseif sls(k)<59
                        Hk(s:f,8) = 1;
                        Hk(s:f,9) = 0;
                    else
                        Hk(s:f,8) = 0;
                        Hk(s:f,9) = 1;
                    end
                elseif sit==10
                    if sls(k)<33
                        Hk(s:f,8) = 0;
                        Hk(s:f,9) = 0;
                    elseif sls(k)<89
                        Hk(s:f,8) = 1;
                        Hk(s:f,9) = 0;
                    else
                        Hk(s:f,8) = 0;
                        Hk(s:f,9) = 1;
                    end
                elseif sit==11
                    if sls(k)<59
                        Hk(s:f,8) = 0;
                        Hk(s:f,9) = 0;
                    elseif sls(k)<89
                        Hk(s:f,8) = 1;
                        Hk(s:f,9) = 0;
                    else
                        Hk(s:f,8) = 0;
                        Hk(s:f,9) = 1;
                    end
                elseif sit==12
                    if sls(k)<33
                        Hk(s:f,8) = 0;
                        Hk(s:f,9) = 0;
                        Hk(s:f,10) = 0;
                    elseif sls(k)<59
                        Hk(s:f,8) = 1;
                        Hk(s:f,9) = 0;
                        Hk(s:f,10) = 0;
                    elseif sls(k)<89
                        Hk(s:f,8) = 0;
                        Hk(s:f,9) = 1;
                        Hk(s:f,10) = 0;
                    elseif sls(k)<106
                        Hk(s:f,8) = 0;
                        Hk(s:f,9) = 0;
                        Hk(s:f,10) = 1;
                    end
                end
                
        end
        
        if sls(k)<33
            Zk(s,1) = i_free(meas((4*k - 3),6),meas((4*k - 2),6),0);
            Zk(f,1) = i_free(meas((4*k - 1),6),meas((4*k    ),6),0);
        elseif sls(k)<59
            Zk(s,1) = i_free(meas((4*k - 3),6),meas((4*k - 2),6),1);
            Zk(f,1) = i_free(meas((4*k - 1),6),meas((4*k    ),6),1);
        elseif sls(k)<89
            Zk(s,1) = i_free(meas((4*k - 3),6),meas((4*k - 2),6),2);
            Zk(f,1) = i_free(meas((4*k - 1),6),meas((4*k    ),6),2);
        elseif sls(k)<106
            Zk(s,1) = i_free(meas((4*k - 3),6),meas((4*k - 2),6),3);
            Zk(f,1) = i_free(meas((4*k - 1),6),meas((4*k    ),6),3);   
        end
        
        p1c = (meas((4*k - 3),7));
        p2c = (meas((4*k - 2),7));
        l1c = (meas((4*k - 1),7));
        l2c = (meas((4*k    ),7));
        if sls(k)<33
            pc = i_free(p1c,p2c,0);
            lc = i_free(l1c,l2c,0);
        elseif sls(k)<59
            pc = i_free(p1c,p2c,1);
            lc = i_free(l1c,l2c,1);
        elseif sls(k)<89
            pc = i_free(p1c,p2c,2);
            lc = i_free(l1c,l2c,2);
        elseif sls(k)<106
            pc = i_free(p1c,p2c,3);
            lc = i_free(l1c,l2c,3);
        end
        
        switch options.TroGrad
            case 0
                if sit==1
                    Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
                    Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));
                elseif sit==2
                    Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
                    Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));
                elseif sit==3 || sit==4 || sit==5
                    if sls(k)<33
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));
                    else
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));
                    end
                elseif sit==6 || sit==7
                    if sls(k)<59
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));
                    else
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));
                    end
                elseif sit==8 || sit==9
                    if sls(k)<33
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));
                    elseif sls(k)<59
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));
                    else
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jn*x1k(bp+k));
                    end
                elseif sit==10
                    if sls(k)<33
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));
                    elseif sls(k)<89
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));
                    else
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jn*x1k(bp+k));
                    end
                elseif sit==11
                    if sls(k)<59
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));
                    elseif sls(k)<89
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));
                    else
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jn*x1k(bp+k));
                    end
                elseif sit == 12
                    if sls(k)<33
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));
                    elseif sls(k)<59
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));
                    elseif sls(k)<89
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Je*x1k(7));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Je*x1k(7)) + (Jn*x1k(bp+k));
                    elseif sls(k)<106
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jc*x1k(8));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jc*x1k(8)) + (Jn*x1k(bp+k));
                    end
                end
            case 1
                if sit==1
                    Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7));
                    Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));
                elseif sit==2
                    Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7));
                    Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));
                elseif sit==3 || sit==4 || sit==5
                    if sls(k)<33
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));
                    else
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));
                    end
                elseif sit==6 || sit==7
                    if sls(k)<59
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));
                    else
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));
                    end
                elseif sit==8 || sit==9
                    if sls(k)<33
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));
                    elseif sls(k)<59
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));
                    else
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jn*x1k(bp+k));
                    end
                elseif sit==10
                    if sls(k)<33
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));
                    elseif sls(k)<89
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));
                    else
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jn*x1k(bp+k));
                    end
                elseif sit==11
                    if sls(k)<59
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));
                    elseif sls(k)<89
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));
                    else
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jn*x1k(bp+k));
                    end
                elseif sit == 12
                    if sls(k)<33
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));
                    elseif sls(k)<59
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));
                    elseif sls(k)<89
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Je*x1k(9));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Je*x1k(9)) + (Jn*x1k(bp+k));
                    elseif sls(k)<106
                        Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jc*x1k(10));
                        Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jc*x1k(10)) + (Jn*x1k(bp+k));
                    end
                end
        end
        
        if strcmp(options.WeMethod,'Elevation Dependent')
            if sls(k)<33
                f1 = freq(sls(k),1); f2 = freq(sls(k),2);
                ab1 = ((f1^2)/(f1^2 - f2^2))^2;
                ab2 = ((f2^2)/(f1^2 - f2^2))^2;
                elv = meas((4*k - 1),26);
                stdcode = (options.CodeStd^2 )*(ab1+ab2);
                stdphas = (options.PhaseStd^2)*(ab1+ab2);
                if elv<30
                    Rk(s,s) = stdcode/(sind(elv));
                    Rk(f,f) = stdphas/(sind(elv));
                else
                    Rk(s,s) = stdcode/(sind(elv));
                    Rk(f,f) = stdphas/(sind(elv));
                end
            elseif sls(k)<59
                f1 = freq(sls(k),1); f2 = freq(sls(k),2);
                ab1 = ((f1^2)/(f1^2 - f2^2))^2;
                ab2 = ((f2^2)/(f1^2 - f2^2))^2;
                elv = meas((4*k - 1),26);
                stdcode = (options.CodeStd^2 )*(ab1+ab2);
                stdphas = (options.PhaseStd^2)*(ab1+ab2);
                if sit==2 || sit==6 || sit==7 || sit==11
                    nr = 1;
                else
                    nr = 4;
                end
                if elv<30
                    Rk(s,s) = stdcode/(sind(elv))*nr;
                    Rk(f,f) = stdphas/(sind(elv));
                else
                    Rk(s,s) = stdcode/(sind(elv))*nr;
                    Rk(f,f) = stdphas/(sind(elv));
                end
            elseif sls(k)<89
                f1 = freq(sls(k),1); f2 = freq(sls(k),2);
                ab1 = ((f1^2)/(f1^2 - f2^2))^2;
                ab2 = ((f2^2)/(f1^2 - f2^2))^2;
                elv = meas((4*k - 1),26);
                stdcode = (options.CodeStd^2 )*(ab1+ab2);
                stdphas = (options.PhaseStd^2)*(ab1+ab2);
                if elv<30
                    Rk(s,s) = stdcode/(sind(elv))*4;
                    Rk(f,f) = stdphas/(sind(elv))*4;
                else
                    Rk(s,s) = stdcode/(sind(elv))*4;
                    Rk(f,f) = stdphas/(sind(elv))*4;
                end
            elseif sls(k)<106
                f1 = freq(sls(k),1); f2 = freq(sls(k),2);
                ab1 = ((f1^2)/(f1^2 - f2^2))^2;
                ab2 = ((f2^2)/(f1^2 - f2^2))^2;
                elv = meas((4*k - 1),26);
                stdcode = (options.CodeStd^2 )*(ab1+ab2);
                stdphas = (options.PhaseStd^2)*(ab1+ab2);
                if elv<30
                    Rk(s,s) = stdcode/(sind(elv))*4;
                    Rk(f,f) = stdphas/(sind(elv))*4;
                else
                    Rk(s,s) = stdcode/(sind(elv))*4;
                    Rk(f,f) = stdphas/(sind(elv))*4;
                end
            end
        end
    end
    
    Ik  = Zk - Ck;
    
    if rcond(Rk)~= 0
        N  = Hk'*(Rk\Hk);
        tx = (Hk'*(Rk\Ik));
    else
        N  = Hk'*(pinv(Rk)*Hk);
        tx = (Hk'*(pinv(Rk)*Ik));
    end
    if rcond(N)~=0
        dx = N\tx;
    else
        dx = pinv(N)*tx;
    end
    
    if norm(dx(1:3,1))>0.001
        x1k = x1k + dx;
    else
        xk = x1k + dx;
        pk  = pinv(N);
        kof = pinv(N);
        kof = kof(1:5,1:5); 
        break
    end
end
end

