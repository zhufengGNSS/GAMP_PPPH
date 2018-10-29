function [] = out_write(xs,Option,model,filename,data)


year = data.inf.time.first(1);
doy  = data.inf.time.doy;

dat = xs';

td = model(1,31);
tw = dat(:,5);
tt = td + tw;

name = {'Year','DOY','SOD','(X)','(Y)','(Z)','(DT)','(TH)','(TW)','(TT)'};
format1 = '%4s %3s %5s %13s %13s %13s %10s %7s %7s %7s';
format2 = '%4d %3d %5d %13.3f %13.3f %13.3f %10.3f %7.3f %7.3f %7.3f';

if Option.TroGrad==1
    nn = size(name,2);
    name{nn+1} = '(TGN)';
    name{nn+2} = '(TGE)';
    
    format1 = [format1,' %7s %7s'];
    format2 = [format2,' %7.4f %7.4f'];
end

if Option.system.glo==1
    nn = size(name,2);
    name{nn+1} = '(SDR)';
    
    format1 = [format1,' %7s'];
    format2 = [format2,' %7.3f'];
end

if Option.system.gal==1
    nn = size(name,2);
    name{nn+1} = '(SDE)';
    
    format1 = [format1,' %7s'];
    format2 = [format2,' %7.3f'];
end

if Option.system.bds==1
    nn = size(name,2);
    name{nn+1} = '(SDC)';
    
    format1 = [format1,' %7s'];
    format2 = [format2,' %7.3f'];
end

fid = fopen(filename,'wt');

fprintf(fid,[format1,'\n'],name{:});
sm = Option.system.glo + Option.system.gal + Option.system.bds;

for i=1:size(dat,1)
    if Option.TroGrad==1
        if sm==0
            fprintf(fid,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6:7));
        elseif sm==1
            fprintf(fid,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6:8));
        elseif sm==2
            fprintf(fid,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6:9));
        elseif sm==3
            fprintf(fid,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6:10));
        end     
    else
        if sm==0
            fprintf(fid,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1));
        elseif sm==1
            fprintf(fid,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6));
        elseif sm==2
            fprintf(fid,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6:7));
        elseif sm==3
            fprintf(fid,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6:8));
        end
    end
end
end

