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

classdef PDE1dImpl < handle

  properties
     xmesh, tspan;
		 pdeFunc,icFunc,bcFunc;
     numDepVars;
     mCoord;
   end
  
  methods
    
    function obj=PDE1dImpl(m, pdeFunc,icFunc,bcFunc,x,t,opts)
    if(~isscalar(m) || ~isnumeric(m) || (m ~= 0 && m ~=1 && m ~=2))
      error(['m must equal 0, 1, or 2 for Cartesian, cylindrical' ...
     ' or spherical coordinate systems, respectively.']);
    end
    obj.mCoord = m;
    obj.pdeFunc = pdeFunc;
		obj.icFunc = icFunc;
		obj.bcFunc = bcFunc;
		obj.xmesh = x;
    obj.numNodes = length(x);
		obj.tspan = t;
    if ~isa(opts, 'PDEOptions')
      error('Argument seven must be a PDEOPtions instance.');
    end
    obj.pdeOpts = opts;
    obj.intRule = GaussianIntegrationRule(GaussianIntegrationRule.Curve, ...
      opts.numIntegrationPoints);
    ic1 = icFunc(x(1));
    obj.numDepVars = length(ic1);
    obj.numElemNodes = 2;
    obj.numFEMEqns = obj.numNodes*obj.numDepVars;
    obj.totalNumEqns = obj.numFEMEqns;
    obj.hasODE = opts.hasODE;
    obj.addLagMultVector = opts.addLagMultVector;
    obj.useDiagMassMat = opts.useDiagMassMat;
    obj.vectorized = opts.vectorized;
    obj.numIntegrationPoints = opts.numIntegrationPoints;
    obj.icDiagnostics = opts.icDiagnostics;
    obj.eqnDiagnostics = opts.eqnDiagnostics;
    obj.analyticalJacobian = opts.analyticalJacobian;

    obj.y0FEM = zeros(obj.numDepVars, obj.numNodes);
    for i=1:obj.numNodes
      obj.y0FEM(:,i) = icFunc(x(i));
    end
    obj.y0 = obj.y0FEM(:);
    
    ne = obj.numNodes-1;
    numXpts = ne*obj.numIntegrationPoints;
    % calc xpts for the mesh
    gPts = obj.intRule.points;
    N = obj.shapeLine2(gPts);
    x = obj.xmesh;
    ip = 1:obj.numIntegrationPoints:numXpts;
    for i=1:obj.numIntegrationPoints
      obj.xPts(ip) = x(1:end-1)*N(1,i) + x(2:end)*N(2,i);
      ip = ip + 1;
    end
    
    obj.isOctave = exist('OCTAVE_VERSION', 'builtin');
    
    end
    
    function setODE(self, odeFunc, odeICFunc, odeMesh)
      self.odeImpl = ODEImpl(self.xmesh, odeFunc, odeICFunc, ...
        odeMesh);
      self.y0 = [self.y0FEM(:); self.odeImpl.y0Ode(:)];
      self.totalNumEqns = size(self.y0,1);
    end
		
    function [outTimes,u,varargout]=solveTransient(self, odeOpts)
      
      reltol=odeOpts.RelTol;
      abstol=odeOpts.AbsTol;
      % flag dirichlet contraints
      xl = self.xmesh(1);
      xr = self.xmesh(end);
      ul = self.y0FEM(:,1);
      ur = self.y0FEM(:,end);
      t0 = self.tspan(1);
      y0Ode = self.odeVFromSysVec(self.y0);
      if self.hasODE
        yp0Ode = zeros(self.odeImpl.numODEDOFs,1);
      else
        yp0Ode = [];
      end
      [pl,ql,pr,qr] = callVarargFunc(self.bcFunc,...
        {xl, ul,xr,ur,t0,y0Ode,yp0Ode});
      self.dirConsFlagsLeft = ql~=0;
      self.dirConsFlagsRight = qr~=0;
      %fprintf('numDepVars=%d\n', obj.numDepVars);
      %prtShortVec(obj.y0, 'y0');
      
      massFunc = @(time, u, up) self.calcMassMat(time, u, up);
      rhsFunc = @(time, u, up) self.calcRHS(time, u, up);
      %opts=odeset('reltol', reltol, 'abstol', abstol);
      %opts=odeset(opts, 'Stats','on');
      opts=odeset(odeOpts);
      icdiag = self.icDiagnostics;
      useDecic=true;
      useOde15i=true;
      if(useOde15i)
        useDecic = true;
      end
      if(useDecic)
        %icf = @(t,y,yp) massFunc(t,y,yp)*yp-rhsFunc(t,y,yp);
        depVarClassType=self.y0;
        icf = @(t,y,yp) self.calcResidual(t,y,yp,depVarClassType);
        t0 = self.tspan(1);
        n=length(self.y0);
        %fixed_y0 = ones(n,1);
        %fixed_y0 = [ones(n-1,1);0];
        %fixed_y0 = [0; ones(n-1,1)];
        fixed_y0 = zeros(n,1);
        %fixed_y0 = [ones(self.numFEMEqns,1);zeros(self.numODE,1)];
        yp0 = zeros(n,1);
        %icf(0,y0,yp0)'
        fixed_yp0 = zeros(n,1);
        icopts.AbsTol=abstol;
        icopts.RelTol=reltol;
        icopts.icdiagnostics=icdiag;
        y0Save = self.y0;
        if(0)
          [y0,yp0]=decic(icf, t0,self.y0,fixed_y0,yp0,fixed_yp0);
        else
          [y0,yp0]=decicShampine(icf, t0,self.y0,fixed_y0,yp0,fixed_yp0,icopts);
        end
        if(icdiag)
          fprintf('max change, y0=%g\n', max(abs(y0-y0Save)));
          prtShortVec(y0, 'y0');
          prtShortVec(yp0, 'yp0', 1, '%16.10g');
        end
        opts=odeset(opts,'initialslope', yp0);
        if(icdiag)
          prtShortVec(icf(t0,y0,yp0), 'res');
        end
      end
      
      abstolVec=abstol*ones(self.totalNumEqns,1);
      % use large convergence tolerance on algebraic eqns
      M= self.calcMassMat(t0, y0, yp0);
      algEqns = abs(diag(M))<eps;
      abstolVec(algEqns) = 1e6*abstol;
      %prtShortVec(abstolVec, 'abstolVec');
      opts=odeset(opts, 'abstol', abstolVec);
      if(self.analyticalJacobian)
        %disp('Using MatlabAutoDiff for Jacobian Evaluations');
        jacFunc = @(time, u, up) self.calcAnalJacobianODE(time, u, up);
        opts=odeset(opts, 'jacobian', jacFunc);
      elseif self.pdeOpts.useInternalNumJac
        jacFunc = @(time, u, up) self.calcSparseJacobian(time, u, up);
        opts=odeset(opts, 'jacobian', jacFunc);
      else 
        % octave ode15i doesn't have JPattern option
        jPat = self.calcJacPattern;
        opts=odeset(opts, 'jpattern', {jPat, jPat});
      end
      
      if(useOde15i)
        %sol=ode15i(icf, self.tspan, y0, yp0, opts);
        [outTimes,u]=ode15i(icf, self.tspan, y0, yp0, opts);
      else
        opts=odeset(opts, 'Mass', massFunc);
        [outTimes,u]=ode15s(rhsFunc, self.tspan, y0, opts);
      end
      if self.hasODE
        uOde = u(:,self.numFEMEqns+1:end);
        varargout{1} = uOde;
        u= u(:,1:self.numFEMEqns);
      end
    end
    
    function [u,varargout]=solveStatic(self, odeOpts)
      
      reltol=odeOpts.RelTol;
      abstol=odeOpts.AbsTol;
      tol=1e-12;
      % flag dirichlet contraints
      xl = self.xmesh(1);
      xr = self.xmesh(end);
      ul = self.y0FEM(:,1);
      ur = self.y0FEM(:,end);
      t0 = self.tspan(1);
      y0Ode = self.odeVFromSysVec(self.y0);
      yp0Ode = 0*y0Ode;
      [pl,ql,pr,qr] = callVarargFunc(self.bcFunc, ...
        {xl, ul,xr,ur,t0,y0Ode,yp0Ode});
      self.dirConsFlagsLeft = ql~=0;
      self.dirConsFlagsRight = qr~=0;
      %fprintf('numDepVars=%d\n', obj.numDepVars);
      %prtShortVec(obj.y0, 'y0');
      
      u=self.y0;
      up=zerosLike(self.totalNumEqns,1,u);
      maxIter = 10;
      it=1;
      converged = false;
      debug = 0;
      resFunc = @(u) self.calcResidual(0,u,up,u);
      while(it<maxIter)
        if debug> 1
        self.printSystemVector(u, 'u');
        end
        r=resFunc(u);
        rn = max(abs(r));
        if debug
          fprintf('it=%2d, rn=%g\n', it, rn);
        end
        if(rn < tol)
          converged = true;
          break;
        end
        if self.analyticalJacobian
          J=self.calcAnalJacobianODE(0, u, up);
        else
          J=self.calcJacobian(0, u, up);
        end      
        %prtMat(J,'J');
        du = J\r;
        if debug> 1
        self.printSystemVector(du, 'du');
        end
        u = u - du;
        it = it + 1;
      end
      if ~converged
        warning('pde1d:maxiter', 'Newton iteration did not converge.');
      end
      if(self.hasODE)
        uODE = u(:,end-self.numODE+1:end);
        varargout{1} = uODE;
      end  
      u=self.femU2FromSysVec(u);
    end
    
    function testODECalc(self)
      ne = self.totalNumEqns;
      if 1
        u = linspace(2,10,ne);
        up = linspace(-2,-4,ne);
      else
        u=0:ne;
        up=u;
      end
      v=self.odeVFromSysVec(u);
      vDot=self.odeVFromSysVec(up);
      [F, S,Cv] = calcFEMEqns(self, 0, u, up, v, vDot);
      prtShortVec(F, 'F');
      prtShortVec(Cv, 'Cxd');
      R = Cv - S + F;
      prtShortVec(R, 'R');
      u2 = self.femU2FromSysVec(u);
      up2 = self.femU2FromSysVec(up);
      f2 = self.femU2FromSysVec(F);
      f=self.odeImpl.calcODEResidual(time, u2, up2, f2, v, vDot);
      dFdu=self.odeImpl.calcDOdeDu(time, u2, up2, f2, v, vDot);
    end
    
    function testFuncs(self)
      ne = self.totalNumEqns;
      %u = zeros(ne, 1);
      u=self.y0;
      up = zeros(ne, 1);
      %up=((1:ne)-1)';
      self.printSystemVector(u, 'u', 1, '%16.8e');
      if self.hasODE
        odeIndices =self.numFEMEqns+1:ne;
        v = u(odeIndices);
        vDot = up(odeIndices);
        uFEM = u(1:self.numFEMEqns);
        upFEM = up(1:self.numFEMEqns);
      else
        v = [];
        vDot = [];
        uFEM = u;
        upFEM = up;
      end
      dudtFunc = self.pdeOpts.initialSlope;
      if ~isempty(dudtFunc)
        dudt2FEM = zeros(self.numDepVars, self.numNodes);
        for i=1:self.numNodes
          dudt2FEM(:,i) = dudtFunc(self.xmesh(i));
        end
        upFEM = dudt2FEM(:);
        up(1:self.numFEMEqns) = upFEM;
      end
      if ~isempty(self.pdeOpts.eqnDiagnosticsInitFunc)
        [uFEM, upFEM] = self.pdeOpts.eqnDiagnosticsInitFunc();
        u(1:self.numFEMEqns) = uFEM;
        up(1:self.numFEMEqns) = upFEM;
      end
      self.printSystemVector(up, 'up', 1, '%16.8e');
      analJac = self.analyticalJacobian;
      if(analJac && self.eqnDiagnostics<1)
        autoDiffF = @(u) self.calcResidual(0, u, up, u);
        for i=1:10
          kAD = AutoDiffJacobianAutoDiff(autoDiffF, u);
        end
        size(kAD)
        return
      end
      depVarClassType=u;
      [F, S,Cxd] = calcFEMEqns(self, 0, uFEM, upFEM, v, vDot, depVarClassType);
      self.printSystemVector(S, 'S', 1, '%16.8e');
      self.printSystemVector(F, 'F', 1, '%16.8e');
      self.printSystemVector(Cxd, 'Cxd', 1, '%16.8e');
      if(self.eqnDiagnostics>1)
        K = self.calcJacobian(0, u, up);
        prtMat(K, 'K', 1, '%13g');
        fprintf('K: rank=%d, cond num=%g\n', rank(K), cond(K));
        M = self.calcMassMat(0, u, up);
        prtMat(M, 'M', 1, '%13g');
        Mup = M*up;
        self.printSystemVector(Mup, 'Mup');
        jPattern = self.calcJacPattern;
        prtMat(full(jPattern), 'jPattern', 1, ' %g');
        if(analJac)
          autoDiffF = @(u) self.calcResidual(0, u, up, u);
          kAD = AutoDiffJacobianAutoDiff(autoDiffF, u);
          prtMat(kAD, 'kAD', 1, '%13g');
        end
      end
      R = self.calcResidual(0, u, up, depVarClassType);
      self.printSystemVector(R, 'R', 1, '%20.10g');
      if self.hasODE
        u2 = self.femU2FromSysVec(u);
        up2 = self.femU2FromSysVec(up);
        f2 = self.femU2FromSysVec(F);        
        [dFdv,dFdvDot]=self.odeImpl.calcDOdeDv(0, u2, up2, f2, v, vDot);
        prtMat(dFdv, 'dFdv', 1, '%10g');
        prtMat(dFdvDot, 'dFdvDot', 1, '%10g');
      end
    end
    
  end % methods
	
	methods(Access=private)

    function jac=calcJacobian(self, time, u, up)
      u=u(:);
      up=up(:);
      ne = length(u);
      jac = zeros(ne,ne);
      depVarClassType=u;
      r0 = calcResidual(self, time, u, up, depVarClassType);
      sqrtEps = sqrt(eps);
      for i=1:ne
        usave = u(i);
        h = sqrtEps * max(usave, 1);
        u(i) = u(i) + h;
        rp = calcResidual(self, time, u, up, depVarClassType);      
        jac(:,i) = (rp-r0)/h;     
        u(i) = usave;
      end
    end
    
    function jac=calcJacobianCD(self, time, u)
      jac = zeros(self.numFEMEqns, self.numFEMEqns);
      disp('calc dFdu by central difference');
      depVarClassType=u;
      for i=1:self.numFEMEqns
        usave = u(i);
        delta = 1e-4;
        u(i) = u(i) + delta;
        rp = calcRHS(self, time, u, depVarClassType);
        u(i) = usave;
        u(i) = u(i) - delta;
        rm = calcRHS(self, time, u, depVarClassType);
        jac(:,i) = (rp-rm)/(2*delta);
        u(i) = usave;
      end
    end
    
    function M=calcMassMat(self, time, u, up)
      u=u(:);
      up=up(:);
      ne = length(u);
      M = zeros(ne,ne);
      depVarClassType=up;
      r0 = calcResidual(self, time, u, up, depVarClassType);
      sqrtEps = sqrt(eps);
      for i=1:ne
        upsave = up(i);
        h = sqrtEps * max(upsave, 1);
        up(i) = up(i) + h;
        rp = calcResidual(self, time, u, up, depVarClassType);      
        M(:,i) = (rp-r0)/h;     
        up(i) = upsave;
      end
    end
    
    function R = calcResidual(self, time, u, up, depVarClassType)
      [rhs,Cxd] = calcRHS(self, time, u, up, depVarClassType);
      R = Cxd+rhs;
    end
    
    function [R,Cxd] = calcRHS(self, time, u, up, depVarClassType)
      if self.hasODE
        v = self.odeVFromSysVec(u);
        vDot = self.odeVFromSysVec(up);
        uFEM = self.femUFromSysVec(u);
        upFEM = self.femUFromSysVec(up);
      else
        v = [];
        vDot = [];
        uFEM = u;
        upFEM = up;
      end
      [F, S,Cxd] = calcFEMEqns(self, time, uFEM, upFEM, v, vDot, depVarClassType);
      R = F-S;
      %R = Cxd-R;
      %R=-R;
      if self.hasODE
        u2 = self.femU2FromSysVec(u);
        up2 = self.femU2FromSysVec(up);
        f2 = self.femU2FromSysVec(F);
        f=self.odeImpl.calcODEResidual(time, u2, up2, f2, v, vDot);
        if(self.addLagMultVector)
          dFdu=self.odeImpl.calcDOdeDu(time, u2, up2, f2, v, vDot);
          R = R + dFdu'*v; % add constraint contribution
        end
        R=[R(:);f(:)];
        % f contains the total residual from ODE equations
        % (there is no separate inertia term)
        Cxd=[Cxd; zeros(self.odeImpl.numODEDOFs,1)];
      end
      % add constraints
      [R,Cxd]=self.applyConstraints(R,Cxd, uFEM,v,vDot,time);
      R=-R;
    end 
    
    function [dfdy, dfdyp]=calcAnalJacobianODE(self, time, u, up)
      autoDiffF = @(u) self.calcResidual(time, u, up, u);
      dfdy = AutoDiffJacobianAutoDiff(autoDiffF, u);
      autoDiffF = @(up) self.calcResidual(time, u, up, up);
      dfdyp = AutoDiffJacobianAutoDiff(autoDiffF, up);
    end
    
    function [R,Cxd]=applyConstraints(self, R, Cxd, u,v, vDot,time)
      u2=self.femU2FromSysVec(u);
      x = self.xmesh;
      xl = x(1);
      xr = x(end);
      %[pl,ql,pr,qr] = self.bcFunc(x(1), u2(:,1),x(end),u2(:,end),time);
      [pl,ql,pr,qr] = callVarargFunc(self.bcFunc, ...
        {xl, u2(:,1),xr,u2(:,end),time,v, vDot});
      rightDofOff = self.numFEMEqns - self.numDepVars;
      m=self.mCoord;
      sing = m~=0 && xl==0;
      for i=1:self.numDepVars
        if(ql(i) ~= 0)
          qli = ql(i);
          if ~sing
            if (m == 1)
              qli = qli/xl;
            elseif(m==2)
              qli = qli/(xl*xl);
            end
          end
          R(i) =  R(i) - pl(i)/qli;
        else
          % apply dirichlet constraint
          R(i) = -pl(i);
          Cxd(i) = 0;
        end
        if(qr(i) ~= 0)
          qri = qr(i);
          if (m == 1)
            qri = qri/xr;
          elseif(m==2)
            qri = qri/(xr*xr);
          end
          R(i+rightDofOff) = R(i+rightDofOff) + pr(i)/qri;
        else
          % apply dirichlet constraint using lagrange multiplier
         R(i+rightDofOff) = -pr(i);
         Cxd(i+rightDofOff) = 0;
       end
     end
     %R = -F + S;
     R = -R;
    end
    
    function [F, S, Cv] = calcFEMEqns(self, t, u,up, v, vDot, depVarClassType)
      ndv=self.numDepVars;
      nn  = self.numNodes;
      Cv = zerosLike(ndv, nn, depVarClassType);
      F = zerosLike(ndv, nn, depVarClassType);
      S = zerosLike(ndv, nn, depVarClassType);
      
      gPts = self.intRule.points;
      gWts = self.intRule.wts;
      numIntPts = length(gPts);
      N = self.shapeLine2(gPts);
      dN = self.dShapeLine2(gPts);
      x = self.xmesh(:)';  
      u2 = self.femU2FromSysVec(u);
      up2 = self.femU2FromSysVec(up);
			numXpts = size(self.xPts,2);
      jac = diff(x)/2;
      nen=self.numElemNodes;
      NN = zeros(nen, nen, numIntPts);
      u21 = u2(:,1:end-1);
      u22 = u2(:,2:end);
      uPts = zerosLike(ndv, numXpts,  depVarClassType);
      duPts = zerosLike(ndv, numXpts, depVarClassType);
      ip = 1:numIntPts:numXpts;
      for i=1:numIntPts
        NN(:,:,i) = N(:,i)*N(:,i)';
        % assume two node elements
        uPts(:,ip) = u21*N(1,i) + u22*N(2,i);
        duPts(:,ip) = u21*dN(1,i)./jac + u22*dN(2,i)./jac;
        ip = ip + 1;
      end
      if(self.vectorized)
        % call pde func for all points
        [c,f,s] = callVarargFunc(self.pdeFunc, ...
          {self.xPts,t,uPts,duPts,v, vDot});
        isFullCMat=self.checkCCoefficientMatrixSizeVec(c);
        self.checkCoefficientMatrixSizeVec(f, 'f');
        self.checkCoefficientMatrixSizeVec(s, 's');
        %prtMat(s, 's');
      else
        f = zerosLike(ndv, numXpts, depVarClassType);
        s = zerosLike(ndv, numXpts, depVarClassType);
        for i=1:size(self.xPts,2)
          [ci,fi,si] = ...
            callVarargFunc(self.pdeFunc, ...
            {self.xPts(i),t,uPts(:,i),duPts(:,i),v, vDot});
        isFullCMat=self.checkCCoefficientMatrixSize(ci);
        self.checkCoefficientMatrixSize(fi, 'f');
        self.checkCoefficientMatrixSize(si, 's');
        f(:,i) = fi;
        s(:,i) = si;
        if i==1
          if isFullCMat
            c = zerosLike(ndv, ndv, numXpts, depVarClassType);
          else
            c = zerosLike(ndv, numXpts, depVarClassType);
          end
        end
          if isFullCMat
            c(:,:,i) = ci;
          else
            c(:,i) = ci;
          end
        end
      end
      if(all(c==0))
        error('pde1d:no_parabolic_eqn', ...
          "At least one of the entries in the c-coefficient vector must be non-zero.");
      end
      
      % form global vectors
      m = self.mCoord;
      if(m==1)
        xm = self.xPts;
      elseif(m==2)
        xm = self.xPts.*self.xPts;
      end
      if m> 0
        if ~isFullCMat
          c = c.*xm;
        else
          c = c.*reshape(xm,1,1,[]);
        end
        f = f.*xm;
        s = s.*xm;
      end
      useDiagMM = self.useDiagMassMat;
      e1 = 1:nn-1;
      e2 = 2:nn;
      ip = 1:numIntPts:numXpts;
      for i=1:numIntPts
        dNdx = dN(:,i)./jac;
        fipJac = f(:, ip).*jac*gWts(i);
        F(:,e1) = F(:,e1) + dNdx(1,:).*fipJac;
        F(:,e2) = F(:,e2) + dNdx(2,:).*fipJac;
        sipJac = s(:, ip).*jac*gWts(i);
        S(:,e1) = S(:,e1) + N(1,i)*sipJac;
        S(:,e2) = S(:,e2) + N(2,i)*sipJac;
        if isFullCMat
          cipJacFull = c(:,:,ip).*(reshape(jac,1,1,[])*gWts(i));
        else
          cipJac = c(:, ip).*jac*gWts(i);
        end
        if(useDiagMM)
          Cv(:,e1) = Cv(:,e1) + up2(:,e1).*cipJac/self.numElemNodes;
          Cv(:,e2) = Cv(:,e2) + up2(:,e2).*cipJac/self.numElemNodes;
        else
          if ~isFullCMat
            Cv(:,e1) = Cv(:,e1) + (NN(1,1,i)*up2(:,e1) + NN(1,2,i)*up2(:,e2)).*cipJac;
            Cv(:,e2) = Cv(:,e2) + (NN(1,2,i)*up2(:,e1) + NN(2,2,i)*up2(:,e2)).*cipJac;
          else
            upie1 = reshape(NN(1,1,i)*up2(:,e1) + NN(1,2,i)*up2(:,e2),ndv,1,[]);
            upie2 = reshape(NN(1,2,i)*up2(:,e1) + NN(2,2,i)*up2(:,e2),ndv,1,[]);
            Cv(:,e1) = Cv(:,e1) + reshape(multiprod(cipJacFull, upie1),ndv,[]);
            Cv(:,e2) = Cv(:,e2) + reshape(multiprod(cipJacFull, upie2),ndv,[]);
          end
        end
      ip = ip + 1;
      end
      F=F(:);
      S=S(:);
      Cv=Cv(:);
    end
    
    function printSystemVector(self, v, name, varargin)
      if ~self.hasODE
        m = reshape(v, self.numDepVars, []);
        prtMat(m, name, varargin{:});
      else
        prtShortVec(v, name, varargin{:});
      end
    end
    
    function checkCoefficientMatrixSize(self, coeffs, coeffsName)
      reqSize=[self.numDepVars, 1];
      sc = size(coeffs);
      % in non-vec mode, allow row or comumn vector
      if length(coeffs) ~= self.numDepVars
        fName = func2str(self.pdeFunc);
        error('pde1d:coeffSize', ['The \"%s\" coefficient matrix ' ...
          'returned from from function \"%s\"\nhas size %d x %d but the ' ...
          'expected size is %d x %d.'], coeffsName, fName, ...
          sc(1), sc(2), reqSize(1), reqSize(2));
      end
    end
    
    function isFullCMat=checkCCoefficientMatrixSize(self, c)
      ndv = self.numDepVars;
      sc = size(c);
      if ~ismatrix(c)
        fName = func2str(self.pdeFunc);
        emsg=sprintf(['The \"c\" coefficient matrix ' ...
          'returned from from function \"%s\"\nhas ' ...
          'an incorrect number of dimensions'], fName);
        error('pde1d:coeffSize', emsg);
      end
      if isvector(c) && length(c)==ndv
        isFullCMat=false;
      elseif all(sc == [ndv, ndv])
        isFullCMat=true;
      else
        fName = func2str(self.pdeFunc);
        emsg=sprintf(['The \"c\" coefficient matrix ' ...
          'returned from from function \"%s\"\nhas '], fName);
        emsg=sprintf([emsg 'size %d x %d but the ' ...
          'expected size is %d x 1.'], sc(1), sc(2), ndv);
        error('pde1d:coeffSize', emsg);
      end
    end

    function checkCoefficientMatrixSizeVec(self, coeffs, coeffsName)
      numXpts = size(self.xPts,2);
      reqSize=[self.numDepVars, numXpts];
      sc = size(coeffs);
      if any(sc ~= reqSize)
        fName = func2str(self.pdeFunc);
        error('pde1d:coeffSize', ['The \"%s\" coefficient matrix ' ...
          'returned from from function \"%s\"\nhas size %d x %d but the ' ...
          'expected size is %d x %d.'], coeffsName, fName, ...
          sc(1), sc(2), reqSize(1), reqSize(2));
      end
    end
    
    function isFullCMat=checkCCoefficientMatrixSizeVec(self, coeffs)
      numXpts = size(self.xPts,2);
      ndv = self.numDepVars;
      sc = size(coeffs);
      dimsc = ndims(coeffs);
      if dimsc==2 && all(sc == [ndv, numXpts])
        isFullCMat=false;
      elseif dimsc==3 && all(sc == [ndv, ndv, numXpts])
        isFullCMat=true;
      else
        fName = func2str(self.pdeFunc);
        emsg=sprintf(['The \"c\" coefficient matrix ' ...
          'returned from from function \"%s\"\nhas '], fName);
        if dimsc == 2
          emsg=sprintf([emsg 'size %d x %d but the ' ...
            'expected size is %d x %d.'], sc(1), sc(2), ndv, numXpts);
        elseif dimsc==3
          emsg=sprintf([emsg 'size %d x %d x %d but the ' ...
            'expected size is %d x %d x %d.'], sc(1), sc(2), sc(3), ndv, ndv, numXpts);
        else
          emsg=[emsg 'an incorrect number of dimensions.'];
        end         
        error('pde1d:coeffSize', emsg);
      end
    end
    
    function J=calcJacPattern(self)
      if ~isempty(self.jacobianPattern)
        J = self.jacobianPattern;
        return;
      end
      ndv = self.numDepVars;
      nn = self.numNodes;
      % jacobian is block-tridiagonal with the block size equal numDepVars
      J = kron( spdiags(ones(nn,3),[-1 0 1],nn,nn), ones(ndv,ndv));
      
      % odes may couple with any other vars
      if self.hasODE
        numODE = self.odeImpl.numODEDOFs;
        J21 = ones(numODE, self.numFEMEqns);
        J22 = ones(numODE, numODE);
        J = [J   J21'
          J21 J22];
      end
      self.jacobianPattern = J;
    end
    
    function [dfdy, dfdyp]=calcSparseJacobian(self, time, u, up)
      jPattern = self.calcJacPattern;
      if isempty(self.jacobianGroupList)
        self.jacobianGroupList=findJacobianGroups(jPattern);
      end
      f = @(time, u) self.calcResidual(time, u, up, u);
      dfdy=numericalJacobian(f, jPattern, u, time, self.jacobianGroupList);
      f = @(time, up) self.calcResidual(time, u, up, up);
      dfdyp=numericalJacobian(f, jPattern, up, time, self.jacobianGroupList);
    end
    
    function v=femUFromSysVec(self,v)
      v=v(:);
      v = v(1:self.numFEMEqns);
    end
    
    function v2=femU2FromSysVec(self,v)
      v2 = reshape(v(1:self.numFEMEqns), self.numDepVars, []);
    end
    
    function vode=odeVFromSysVec(self,v)
      v=v(:);
      vode = v(self.numFEMEqns+1:end);
    end
    
	end % methods
  
  methods(Access=private, Static)
    
    function shp = shapeLine2( r )
      shp = [(1-r)/2 (1+r)/2]';
    end
    
    function ds = dShapeLine2( r )
      ds = [-.5 .5]'*ones(1,length(r));
    end
    
  end % methods
  
  properties(Access=private)
    isOctave;
    intRule;
    numNodes, numElemNodes;
    numFEMEqns, totalNumEqns;
    numIntegrationPoints;
    hasODE, odeImpl;
    dirConsFlagsLeft, dirConsFlagsRight;
    y0FEM, y0;
    addLagMultVector;
    useDiagMassMat;
    vectorized;
    icDiagnostics, eqnDiagnostics, analyticalJacobian;
    pdeOpts;
    jacobianPattern, jacobianGroupList;
    % temporary arrays for vectorized mode
    xPts;
  end
  
end % classdef
