function  [Sx,Sy,F] = my_ampls(beam,pairs)
%Calculate maximum amplitude for each grid point
%Using Root mean square (L2)
%--------------------------------------------------------------------------

%sqr_beam=(beam.^2)'; 
%ampls=sqrt(sqr_beam./length(sqr_beam)); 
%Find maximum 
[F,ind]=max(max((beam'))); 
%find maximum
%use later in color-plot
%define backazimuth and slowness vector
Sx=pairs(ind,1); Sy=pairs(ind,2);
end
