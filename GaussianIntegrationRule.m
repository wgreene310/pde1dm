%  
%   This library is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation; either version 3
%   of the License, or (at your option) any later version.
%  
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   General Public License for more details.
%  
%   You should have received a copy of the GNU General Public License
%   along with this library; if not, visit
%   http://www.gnu.org/licenses/gpl.html or write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
% 
% Copyright (C) 2016-2020 William H. Greene

classdef GaussianIntegrationRule
  
  properties
    points, wts;
  end
  
  properties(Constant)
      Curve=1; Triangle=2; Tetrahedron=3;
  end
  
  methods
    function obj = GaussianIntegrationRule(regionType, numPoints)
      if(regionType == GaussianIntegrationRule.Curve)
        obj = setCurvePoints(obj, numPoints);
      elseif(regionType == GaussianIntegrationRule.Triangle)
        obj=setTriPoints(obj, numPoints);
      elseif(regionType == GaussianIntegrationRule.Tetrahedron)
        obj=setTetPoints(obj, numPoints);
      else
        disp('Illegal regionType requested.');
      end
    end
  end
  
  methods(Access=private)
    
    function self=setCurvePoints(self, numPoints)
      if(numPoints==1)
        self.points = 0;
        self.wts = 2;
      elseif(numPoints==2)
        pt = sqrt(3)/3;
        self.points = [-pt pt]';
        self.wts = [1 1]';
      elseif(numPoints==3)
        pt = sqrt(3/5);
        self.points = [-pt 0 pt]';
        self.wts = [5 8 5]'./9;
      elseif(numPoints==4)
        % static const double c4p0 = 2.0*sqrt(6./5.);
        % static const double c4p1 = sqrt((3. - c4p0)/7.);
        % static const double c4p2 = sqrt((3. + c4p0)/7.);
        % static const double c4wt1 = (18. + sqrt(30.))/36.;
        % static const double c4wt2 = (18. - sqrt(30.))/36.;
        
        p0 = 2*sqrt(6/5);
        p1 = sqrt((3 - p0)/7);
        p2 = sqrt((3 + p0)/7);
        self.points = [ -p2 -p1 p1 p2]';
        wt1 = (18 + sqrt(30))/36;
        wt2 = (18 - sqrt(30))/36;
        self.wts = [wt2 wt1 wt1 wt2]';
      elseif(numPoints==5)
        c5c1 = 2*sqrt(10/7); c5c2 = 13*sqrt(70);
        c5p1 = 1/3*sqrt(5 - c5c1); c5p2 = 1/3*sqrt(5 + c5c1);
        c5w1 = (322+c5c2)/900; c5w2 = (322-c5c2)/900;
        self.points = [0 -c5p1 -c5p2 c5p1 c5p2]';
        self.wts = [128/225 c5w1 c5w2 c5w1 c5w2]';
      elseif(numPoints==8)
        tab=[
          0.3626837833783620	-0.1834346424956498
          0.3626837833783620	0.1834346424956498
          0.3137066458778873	-0.5255324099163290
          0.3137066458778873	0.5255324099163290
          0.2223810344533745	-0.7966664774136267
          0.2223810344533745	0.7966664774136267
          0.1012285362903763	-0.96028985649753638
          0.1012285362903763	0.9602898564975363];
        self.points = tab(:,2)'; self.wts = tab(:,1)';
      elseif(numPoints==12)
        tab=[
          0.2491470458134028	-0.1252334085114689
          0.2491470458134028	0.1252334085114689
          0.2334925365383548	-0.3678314989981802
          0.2334925365383548	0.3678314989981802
          0.2031674267230659	-0.5873179542866175
          0.2031674267230659	0.5873179542866175
          0.1600783285433462	-0.7699026741943047
          0.1600783285433462	0.7699026741943047
          0.1069393259953184	-0.9041172563704749
          0.1069393259953184	0.9041172563704749
          0.0471753363865118	-0.9815606342467192
          0.0471753363865118	0.9815606342467192];
        self.points = tab(:,2)'; self.wts = tab(:,1)';
      elseif(numPoints==18)
        tab=[
          0.1691423829631436	-0.0847750130417353
          0.1691423829631436	0.0847750130417353
          0.1642764837458327	-0.2518862256915055
          0.1642764837458327	0.2518862256915055
          0.1546846751262652	-0.4117511614628426
          0.1546846751262652	0.4117511614628426
          0.1406429146706507	-0.5597708310739475
          0.1406429146706507	0.5597708310739475
          0.1225552067114785	-0.6916870430603532
          0.1225552067114785	0.6916870430603532
          0.1009420441062872	-0.8037049589725231
          0.1009420441062872	0.8037049589725231
          0.0764257302548891	-0.8926024664975557
          0.0764257302548891	0.8926024664975557
          0.0497145488949698	-0.9558239495713977
          0.0497145488949698	0.9558239495713977
          0.0216160135264833	-0.9915651684209309
          0.0216160135264833	0.9915651684209309];
          self.points = tab(:,2)'; self.wts = tab(:,1)';
      else
        error('Illegal number of points in setCurvePoints');
      end
    end
    
    function self=setTriPoints(self, numPoints)
      jacFac = .5;
      if(numPoints==1)
        self.points=[1/3,1/3]';
        self.wts=jacFac;
      elseif(numPoints==3)
        %self.points=[.5,.5; .5,0; 0 .5]';
        self.points=[2/3 1/6; 1/6 2/3; 1/6 1/6]';
        self.wts=[1/3, 1/3, 1/3]*jacFac;
      elseif(numPoints==6)
        g16p = (8.-sqrt(10.) + sqrt(38.-44.*sqrt(2./5.)))/18.;
        g26p = (8.-sqrt(10.) - sqrt(38.-44.*sqrt(2./5.)))/18.;
        w16p = (620 + sqrt(213125.-53320.*sqrt(10.)))/3720.;
        w26p = (620 - sqrt(213125.-53320.*sqrt(10.)))/3720.;
        self.points=[
       1.-2.*g16p, g16p
       g16p, 1.-2.*g16p
       g16p, g16p
       1.-2.*g26p, g26p
       g26p, 1.-2.*g26p
       g26p, g26p]';
       self.wts=[w16p w16p w16p w26p w26p w26p]*jacFac;
      else
        error('Illegal number of points in setTriPoints');
      end
    end
    
    function self=setTetPoints(self, numPoints)
      jacFac = 1/6;
      if(numPoints==1)
        self.points=[.25 .25 .25]';
        self.wts=jacFac*1.0;
      elseif(numPoints==4)
        a = (5+3*sqrt(5))/20;
        b = (5-sqrt(5))/20;
        self.points=[a,b,b; b,a,b; b,b,a; b,b,b]';
        self.wts=jacFac*[.25,.25,.25,.25];
      elseif(numPoints==5)
        o6=1./6.;
        self.points=[.25,.25,.25; .5,o6,o6; o6,.5,o6;
                     o6,o6,.5; o6,o6,o6]';
        n20 = 9./20.;
        self.wts=jacFac*[-4./5. n20 n20 n20 n20];
      elseif(numPoints==15)
        s15=sqrt(15);
         g1=(7-s15)/34;  g2=7/17-g1; g3=(10-2*s15)/40; g4=.5-g3;
         w1=(2665+14*s15)/37800; w2=(2665-14*s15)/37800;
         w3=10/189;
         self.points=[
           1-3*g1, g1, g1;
           g1, 1-3*g1, g1;
           g1, g1, 1-3*g1;
           g1, g1, g1;
           1-3*g2, g2, g2;
           g2, 1-3*g2, g2;
           g2, g2, 1-3*g2;
           g2, g2, g2; %8
           g4, g4, g3;
           g4, g3, g4;
           g4, g3, g3;
           g3, g4, g4;
           g3, g4, g3;
           g3, g3, g4;
           .25,.25,.25; % 15
           ]';
         self.wts=jacFac*[w1,w1,w1,w1,w2,w2,w2,w2,w3,w3,w3,w3,w3,w3,16/135];
      elseif(numPoints==24)
         g1=0.214602871259152029288839219386284991; g11 = 1-3*g1;
         g2=0.040673958534611353115579448956410059; g22 = 1-3*g2;
         g3=0.322337890142275510343994470762492125; g33 = 1-3*g3;
         w1= (85+2*g2*(-319+9*sqrt(5)+624*g2)-638*g3- ...
         24*g2*(-229+472*g2)*g3+96*(13+118*g2*(-1+2*g2))*g3^2+ ...
         9*sqrt(5)*(-1+2*g3))/(13440*(g1-g2)*(g1-g3)*(3-8*g2+ ... 
         8*g1*(-1+2*g2)-8*g3+16*(g1+g2)*g3));
         w2= -(85+2*g1*(-319+9*sqrt(5)+624*g1)-638*g3- ...
         24*g1*(-229+472*g1)*g3+96*(13+118*g1*(-1+2*g1))*g3^2+ ...
         9*sqrt(5)*(-1+2*g3))/(13440*(g1-g2)*(g2-g3)*(3-8*g2+ ...
         8*g1*(-1+2*g2)-8*g3+16*(g1+g2)*g3));
         w3= (85+2*g1*(-319+9*sqrt(5)+624*g1)-638*g2- ...
         24*g1*(-229+472*g1)*g2+96*(13+118*g1*(-1+2*g1))*g2^2+ ...
         9*sqrt(5)*(-1+2*g2))/(13440*(g1-g3)*(g2-g3)*(3-8*g2+ ...
         8*g1*(-1+2*g2)-8*g3+16*(g1+g2)*g3));
         g4=(3-sqrt(5))/12; h4=(5+sqrt(5))/12; p4=(1+sqrt(5))/12;
%         If [i<5, info={{g1,g1,g1,g1},w1};info[[1,i]]= 1-3*g1];
%         If [i>4&&i<9, info={{g2,g2,g2,g2},w2};info[[1,i-4]]=1-3*g2];
%         If [i>8&&i<13,info={{g3,g3,g3,g3},w3};info[[1,i-8]]=1-3*g3];
%         If [i>12,info={{g4,g4,g4,g4},27/560};
%            {j,k}=jk12[[i-12]];info[[1,j]]=h4;info[[1,k]]=p4] ]; 
          self.points=[
            g11 g1 g1; g1 g11 g1; g1 g1 g11; g1 g1 g1;
            g22 g2 g2; g2 g22 g2; g2 g2 g22; g2 g2 g2;
            g33 g3 g3; g3 g33 g3; g3 g3 g33; g3 g3 g3;
            p4 g4 g4; g4 p4 g4; g4 g4 p4;
            h4 g4 g4; g4 h4 g4; g4 g4 h4; 
            g4 p4 h4; p4 h4 g4; h4 g4 p4; 
            g4 h4 p4; p4 g4 h4; h4 p4 g4;
          ]';
          self.wts=jacFac*[w1,w1,w1,w1,w2,w2,w2,w2,w3,w3,w3,w3,repmat(27/560,1,12)]';
      elseif(numPoints==46)
         tab=[
  self.S31(.0396754230703899012650713295393895, .0063971477799023213214514203351730)
  self.S31(.3144878006980963137841605626971483, .0401904480209661724881611584798178)
  self.S31( .1019866930627033000000000000000000, .0243079755047703211748691087719226)
  self.S31( .1842036969491915122759464173489092, .0548588924136974404669241239903914)
  self.S22(.0634362877545398924051412387018983, .0357196122340991824649509689966176)
  self.S211(.0216901620677280048026624826249302, .7199319220394659358894349533527348, ...
            .0071831906978525394094511052198038);
  self.S211(.2044800806367957142413355748727453, .5805771901288092241753981713906204, ...
            .0163721819453191175409381397561191)
];
      self.points = tab(:,1:3)';
      self.wts = jacFac*tab(:,5);
      %sum(tab(:,1:4)')
      else
          error('Illegal number of points in setTetPoints');
      end
    end
end % methods

methods(Access=private, Static)
    
   function xw=S31(a,w)
      b=1-3*a;
      xw=[
        a,a,a,b
        a,a,b,a
        a,b,a,a
        b,a,a,a
      ];
      xw(:,5)=w;
    end
    

    function xw=S22(a,w)
     b=.5-a;
     xw=[
       a,a,b,b
       a,b,a,b
       a,b,b,a
       b,a,a,b
       b,a,b,a
       b,b,a,a
      ];
      xw(:,5)=w;
    end
    function xw=S211(a,b,w)
      c = 1-2*a-b;
     xw=[
        a,a,b,c
        a,a,c,b
        a,b,a,c
        a,b,c,a
        a,c,a,b
        a,c,b,a
        b,a,c,a
        b,a,a,c
        c,a,b,a
        c,a,a,b
        c,b,a,a
        b,c,a,a
      ];
     xw(:,5)=w;
    end
    
  end % methods
  
end


    

