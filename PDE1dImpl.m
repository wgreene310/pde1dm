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
% Copyright (C) 2016-2021 William H. Greene

classdef PDE1dImpl < handle
  
  %% Main implementaton class for PDE solver.

  properties
     xmesh, tspan;
		 pdeFunc,icFunc,bcFunc;
     numDepVars;
     mCoord;
   end
  
  methods
    
    function obj=PDE1dImpl(m, pdeFunc,icFunc,bcFunc,x,t,opts, ...
        odeFunc, odeICFunc, odeMesh)
    if(~isscalar(m) || ~isnumeric(m) || (m ~= 0 && m ~=1 && m ~=2))
      error(['m must equal 0, 1, or 2 for Cartesian, cylindrical' ...
     ' or spherical coordinate systems, respectively.']);
    end
    obj.mCoord = m;
    obj.pdeFunc = pdeFunc;
		obj.icFunc = icFunc;
		obj.bcFunc = bcFunc;
    x=x(:)';
		obj.xmesh = x;
    obj.numNodes = length(x);
    t=t(:)';
		obj.tspan = t;
    if ~isa(opts, 'PDEOptions')
      error('Argument seven must be a PDEOPtions instance.');
    end
    obj.pdeOpts = opts;
    if isempty(opts.numIntegrationPoints)
      opts.numIntegrationPoints = opts.polyOrder+1;
    end
    obj.intRule = GaussianIntegrationRule(GaussianIntegrationRule.Curve, ...
      opts.numIntegrationPoints);
    ic1 = icFunc(x(1));
    obj.numDepVars = length(ic1);
    obj.numElemNodes = opts.polyOrder+1;
    obj.numFEMEqns = obj.numNodes*obj.numDepVars;
    obj.totalNumEqns = obj.numFEMEqns;
    obj.hasODE = opts.hasODE;
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
    
    obj.numElems=(obj.numNodes-1)/(obj.numElemNodes-1);
    numXpts = obj.numElems*obj.numIntegrationPoints;
    % calc xpts for the mesh
    gPts = obj.intRule.points;
    if obj.numElemNodes==2
      obj.sF =  @obj.shapeLine2;
      obj.dSf = @obj.dShapeLine2;
      obj.N = obj.shapeLine2(gPts);
      obj.dN = obj.dShapeLine2(gPts);
    elseif obj.numElemNodes==3
      obj.sF =  @obj.shapeLine3;
      obj.dSf = @obj.dShapeLine3;
      obj.N = obj.shapeLine3(gPts);
      obj.dN = obj.dShapeLine3(gPts);
     elseif obj.numElemNodes==4
      obj.sF =  @obj.shapeLine4;
      obj.dSf = @obj.dShapeLine4;
      obj.N = obj.shapeLine4(gPts);
      obj.dN = obj.dShapeLine4(gPts);
    else
      error('Unsupported number of element nodes.');
    end
    x = obj.xmesh;
    obj.xPts=zeros(1, numXpts);
    obj.NN = zeros(obj.numElemNodes, obj.numElemNodes, obj.numIntegrationPoints);
    ip = 1:obj.numIntegrationPoints:numXpts;
    for i=1:obj.numIntegrationPoints
      obj.NN(:,:,i) = obj.N(:,i)*obj.N(:,i)';
      for j=1:obj.numElemNodes
        jj=obj.numElemNodes-j;
        obj.xPts(ip) = obj.xPts(ip) + x(j:obj.numElemNodes-1:end-jj)*obj.N(j,i);
      end
      ip = ip + 1;
    end
    
    obj.isOctave = exist('OCTAVE_VERSION', 'builtin');
    
    if obj.hasODE
      obj.setODE(odeFunc, odeICFunc, odeMesh);
    end
    
    end
    
		
    function [u,varargout]=solveTransient(self, odeOpts)
      
      reltol=odeOpts.RelTol;
      abstol=odeOpts.AbsTol;
      % flag dirichlet contraints
      xl = self.xmesh(1);
      xr = self.xmesh(end);
      ul = self.y0FEM(:,1);
      ur = self.y0FEM(:,end);
      t0 = self.tspan(1);
      
      if self.hasODE
        y0Ode = self.odeImpl.getVFromSysVec(self.y0);
        yp0Ode = zeros(self.odeImpl.numODEEquations,1);
      else
        y0Ode = [];
        yp0Ode = [];
      end
      [pl,qLeft,pr,qRight] = callVarargFunc(self.bcFunc,...
        {xl, ul,xr,ur,t0,y0Ode,yp0Ode});
      self.dirConsFlagsLeft = qLeft~=0;
      self.dirConsFlagsRight = qRight~=0;
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
        yp0 = self.pdeOpts.initialSlope;
        if isempty(yp0)
          yp0 = zeros(n,1);
        end
        %fixed_y0 = ones(n,1);
        %fixed_y0 = [ones(n-1,1);0];
        %fixed_y0 = [0; ones(n-1,1)];
        fixed_y0 = zeros(n,1); % all y0 free
        if 0
          % this strategy is too restrictive
          % fix all non-BC fem dofs
          fixed_y0(self.numDepVars+1:self.numFEMEqns-self.numDepVars) = 1;
          % free any dofs that are obviously algebraic
          [dfdy,dfdyp]=self.calcSparseJacobian(t0, self.y0, yp0);
          if 0
            for i=1:n
              if all(dfdyp(i,:)==0) && all(dfdyp(:,i)==0)
                fixed_y0(i)=0;
              end
            end
          else
            [Q,R,E] = qr(full(dfdyp));
            adr=abs(diag(R));
            tol=eps*max(adr);
            rr=nnz(adr>tol);
            ii=rr+1:n;
            free_y0=logical(sum(E(:,ii)'));
            fixed_y0(free_y0)=0;
          end
        end
        %icf(0,y0,yp0)'
        fixed_yp0 = zeros(n,1);
        if ~isempty(self.pdeOpts.cicAbsTol)
          icopts.AbsTol= self.pdeOpts.cicAbsTol;
        else
          icopts.AbsTol=abstol/10;
        end
        if ~isempty(self.pdeOpts.cicRelTol)
          icopts.RelTol= self.pdeOpts.cicRelTol;
        else
          icopts.RelTol=reltol/10;
        end
        icopts.Jacobian = @self.calcJacobian;
        
        y0=self.y0;
        y0Save = y0;
        if self.pdeOpts.cicMethod==1
          [y0,yp0]=decic(icf, t0,y0,fixed_y0,yp0,fixed_yp0,icopts);
        elseif self.pdeOpts.cicMethod==0
          icopts.icdiagnostics=icdiag;
          icopts.maxiter=30;
          [y0,yp0]=decicShampine(icf, t0,y0,fixed_y0,yp0,fixed_yp0,icopts);
        end
        if(icdiag)
          prtShortVec(self.xmesh, 'x');
          chgTol=max(icopts.RelTol*abs(self.y0), icopts.AbsTol);
          if 0
          [maxChg, indMax] = max(abs(y0-y0Save));
          fprintf('max change, y0=%g at equation %d\n', maxChg, indMax);
          end
          chgY = y0-y0Save;
          chgI = abs(chgY) > chgTol;
          prtShortVec(chgY(chgI), 'change in y0 values');
          I=1:n;
          prtShortVec(I(chgI), 'changed y0 equations');
          prtShortVec(y0, 'y0');
          prtShortVec(yp0, 'yp0', 1, '%16.8g');
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
        jacFunc = @(time, u, up) self.calcJacobianAutoDiff(time, u, up);
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
        [outTimes,u]=ode15i(icf, self.tspan, y0, yp0, opts);
      else
        opts=odeset(opts,'initialslope', yp0);
        opts=odeset(opts, 'Mass', massFunc);
        [outTimes,u]=ode15s(rhsFunc, self.tspan, y0, opts);
      end
      if self.hasODE
        uOde = u(:,self.numFEMEqns+1:self.numFEMEqns+self.odeImpl.numODEVariables);
        varargout{1} = uOde;
        u= u(:,1:self.numFEMEqns);
      end
    end
    
    function [u,varargout]=solveStatic(self, odeOpts)
      
      reltol=odeOpts.RelTol;
      abstol=odeOpts.AbsTol;
      tol=1e-10;
      % flag dirichlet contraints
      xl = self.xmesh(1);
      xr = self.xmesh(end);
      ul = self.y0FEM(:,1);
      ur = self.y0FEM(:,end);
      t0 = self.tspan(1);
      if 0
      y0Ode = self.odeVFromSysVec(self.y0);
      yp0Ode = 0*y0Ode;
      else
        y0Ode=[];
        yp0Ode=[];
      end
      [pl,ql,pr,qr] = callVarargFunc(self.bcFunc, ...
        {xl, ul,xr,ur,t0,y0Ode,yp0Ode});
      self.dirConsFlagsLeft = ql~=0;
      self.dirConsFlagsRight = qr~=0;
      %fprintf('numDepVars=%d\n', obj.numDepVars);
      %prtShortVec(obj.y0, 'y0');
      
      u=self.y0;
      up=zeros(self.totalNumEqns,1,'like', u);
      maxIter = 10;
      it=1;
      converged = false;
      debug = 1;
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
          J=self.calcJacobianAutoDiff(0, u, up);
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
      [F, S,Cxd, fIntPts] = calcFEMEqns(self, 0, uFEM, upFEM, v, vDot, depVarClassType);
      self.printSystemVector(S, 'S', 1, '%16.8e');
      self.printSystemVector(F, 'F', 1, '%16.8e');
      self.printSystemVector(Cxd, 'Cxd', 1, '%16.8e');
      self.printSystemVector(S-F, 'S-F', 1, '%16.8e');
      self.printSystemVector(S-F+Cxd, 'S-F+Cxd', 1, '%16.8e');
      self.printSystemVector(S-F-Cxd, 'S-F-Cxd', 1, '%16.8e');
      self.printSystemVector(Cxd+S-F, 'Cxd+S-F', 1, '%16.8e');
      if(self.eqnDiagnostics>1)
        K = self.calcJacobian(0, u, up);
        prtMat(K, 'K', 1, '%13g');
        fprintf('K: rank=%d, cond num=%g\n', sprank(K), condest(K));
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
        us=SysVec(u, self.numNodes, self.numDepVars);
        ups=SysVec(up, self.numNodes, self.numDepVars);
        self.odeImpl.testFuncs(us, ups, self.xPts, fIntPts, v, vDot);
      end
    end
    
    function R = calcResidual(self, time, u, up, depVarClassType)
      if nargin<4
        depVarClassType=u;
      end
      [rhs,Cxd] = calcRHS(self, time, u, up, depVarClassType);
      R = Cxd+rhs;
    end
    
    function [R,Cxd] = calcRHS(self, time, u, up, depVarClassType)
      if nargin<4
        depVarClassType=u;
      end
      if self.hasODE
        v = self.odeImpl.getVFromSysVec(u);
        vDot = self.odeImpl.getVFromSysVec(up);
        uFEM = self.femUFromSysVec(u);
        upFEM = self.femUFromSysVec(up);
      else
        v = [];
        vDot = [];
        uFEM = u;
        upFEM = up;
      end
      [F, S,Cxd, fIntPts] = calcFEMEqns(self, time, uFEM, upFEM, v, vDot, depVarClassType);
      R = F-S;
      if self.hasODE
        us=SysVec(u, self.numNodes, self.numDepVars);
        ups=SysVec(up, self.numNodes, self.numDepVars);
        [R,Cxd]=self.odeImpl.updateResiduals(time, us, ups, ...
          self.xPts, fIntPts, R, Cxd);
      end
      % add constraints
      [R,Cxd]=self.applyConstraints(R,Cxd, uFEM,v,vDot,time);
      R=-R;
    end 
    
    function [dfdy, dfdyp]=calcJacobian(self, time, u, up)
      if(self.analyticalJacobian)
        [dfdy, dfdyp]=self.calcJacobianAutoDiff(time, u, up);
      else
        [dfdy, dfdyp]=self.calcSparseJacobian(time, u, up);
      end
    end
    
    function [dfdy, dfdyp]=calcJacobianAutoDiff(self, time, u, up)
      %u=u(:);
      %up=up(:);
      autoDiffF = @(u) self.calcResidual(time, u, up, u);
      dfdy = AutoDiffJacobianAutoDiff(autoDiffF, u);
      autoDiffF = @(up) self.calcResidual(time, u, up, up);
      dfdyp = AutoDiffJacobianAutoDiff(autoDiffF, up);
    end
    
  end % methods
	
  methods(Access=private)
    
    function setODE(self, odeFunc, odeICFunc, odeMesh)
      v0=odeICFunc();
      nV = size(v0,1);
      v0Dot = zeros(nV,1);
      t0 = self.tspan(1);
      yp0 = zeros(self.numFEMEqns, 1);
      [F, S,Cxd,fIntPts] = self.calcFEMEqns(t0, self.y0, yp0, v0, v0Dot, self.y0);
      u2 = self.femU2FromSysVec(self.y0);
      up2 = self.femU2FromSysVec(yp0);
      meshMapper = PDEMeshMapper(self.xmesh, odeMesh, self.numElemNodes, self.sF, self.dSf);
      self.odeImpl = ODEImpl(meshMapper, odeFunc, odeICFunc, ...
        odeMesh, t0, u2, up2, self.xPts, fIntPts, self.pdeOpts.testFunctionDOFMap);
      numODEEqn = self.odeImpl.numODEEquations;
      self.totalNumEqns = self.numFEMEqns + numODEEqn;
      % there can be more ODE equations than ODE variables if
      %  some of the ODE equations are just constraints on the FEM dofs
      self.y0 = [self.y0FEM(:); v0; zeros(numODEEqn-nV,1)];
    end

    function jac=calcJacobianFullFD(self, time, u, up)
      % testing only
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
    
    function jac=calcJacobianFullCD(self, time, u)
      % testing only
      % currently unused
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
    
    function [R,Cxd]=applyConstraints(self, R, Cxd, u,v, vDot,time)
      u2=self.femU2FromSysVec(u);
      x = self.xmesh;
      xl = x(1);
      xr = x(end);
      %[pl,ql,pr,qr] = self.bcFunc(x(1), u2(:,1),x(end),u2(:,end),time);
      [pl,ql,pr,qr] = callVarargFunc(self.bcFunc, ...
        {xl, u2(:,1),xr,u2(:,end),time,v, vDot});
      self.checkBCSize(pl);
      self.checkBCSize(ql);
      self.checkBCSize(pr);
      self.checkBCSize(qr);
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
    
    function [F, S, Cv, f] = calcFEMEqns(self, t, u,up, v, vDot, depVarClassType)
      ndv=self.numDepVars;
      nn  = self.numNodes;
      Cv = zeros(ndv, nn, 'like', depVarClassType);
      F = zeros(ndv, nn, 'like', depVarClassType);
      S = zeros(ndv, nn, 'like', depVarClassType);
      
      gPts = self.intRule.points;
      gWts = self.intRule.wts;
      numIntPts = length(gPts);
      nen=self.numElemNodes;
      N=self.N; %#ok<PROPLC>
      dN=self.dN; %#ok<PROPLC>
      NN=self.NN; %#ok<PROPLC>
      x = self.xmesh(:)';
      u2 = self.femU2FromSysVec(u);
      up2 = self.femU2FromSysVec(up);
			numXpts = size(self.xPts,2);
      nenM1=nen-1;
      i1=1:nenM1:nn-nen+1;
      i2=nen:nenM1:nn;
      jac=(x(i2)-x(i1))/2;
      uPts = zeros(ndv, numXpts,  'like', depVarClassType);
      duPts = zeros(ndv, numXpts, 'like', depVarClassType);
      ip = 1:numIntPts:numXpts-numIntPts+1;
      for i=1:numIntPts
        jj=nenM1;
        for j=1:nen
          u2j=u2(:,j:nenM1:end-jj);
          uPts(:,ip) = uPts(:,ip) + u2j*N(j,i); %#ok<PROPLC>
          duPts(:,ip) = duPts(:,ip) + u2j*dN(j,i); %#ok<PROPLC>
          jj=jj-1;
        end
        duPts(:,ip) = duPts(:,ip)./jac;
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
        f = zeros(ndv, numXpts, 'like', depVarClassType);
        s = zeros(ndv, numXpts, 'like', depVarClassType);
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
            c = zeros(ndv, ndv, numXpts, 'like', depVarClassType);
          else
            c = zeros(ndv, numXpts, 'like', depVarClassType);
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
      ip = 1:numIntPts:numXpts;
      for i=1:numIntPts
        dNdx = dN(:,i)./jac;
        fipJac = f(:, ip).*jac*gWts(i);
        sipJac = s(:, ip).*jac*gWts(i);
        if isFullCMat
          cipJacFull = c(:,:,ip).*(reshape(jac,1,1,[])*gWts(i));
        else
          cipJac = c(:, ip).*jac*gWts(i);
        end
        for j=1:nen
          ej=j:nenM1:nn-nen+j;
          F(:,ej) = F(:,ej) + dNdx(j,:).*fipJac;
          S(:,ej) = S(:,ej) + N(j,i)*sipJac;
          if(useDiagMM)
            Cv(:,ej) = Cv(:,ej) + up2(:,ej).*cipJac/nen;
          else
            for k=1:nen
              ek=k:nenM1:nn-nen+k;
              if ~isFullCMat
                Cv(:,ej) = Cv(:,ej) + NN(j,k,i)*up2(:,ek).*cipJac;
              else
                upiek = reshape(NN(j,k,i)*up2(:,ek) ,ndv,1,[]);
                Cv(:,ej) = Cv(:,ej) + reshape(multiprod(cipJacFull, upiek),ndv,[]);
              end
            end
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
      po = self.pdeOpts.polyOrder;
      inds=-po:po;
      J = kron( spdiags(ones(nn,length(inds)),inds,nn,nn), ones(ndv,ndv));
      
      % odes may couple with any other vars
      if self.hasODE
        numODE = self.odeImpl.numODEEquations;
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
    
    function checkBCSize(self, bc)
      sb = size(bc);
      if any(sb ~= [self.numDepVars 1])
        bcName = inputname(2);
        bcErrMsg = 'Size of "%s" must be %d X 1';
        error('pde1d:bcSize', bcErrMsg, bcName, sb(1));
      end
    end
    
	end % methods
  
  methods(Access=public, Static)
    
    function shp = shapeLine2( r )
      shp = [(1-r)/2 (1+r)/2]';
    end
    
    function ds = dShapeLine2( r )
      ds = [-.5 .5]'*ones(1,length(r));
    end
    
    function shp = shapeLine3( r )
      shp = [-r.*(1-r)/2 1-r.^2 r.*(1+r)/2]';
    end
    
    function ds = dShapeLine3( r )
      r=r(:)';
      ds = [-.5+r; -2*r; .5+r];
    end
    
    function shp = shapeLine4( r )
      r=r(:)';
      rp1=r+1;
      rm1=r-1;
      rp13=r+1/3;
      rm13=r-1/3;
      shp = [-9/16*(rp13).*(rm13).*(rm1)
        27/16*(rp1).*(rm13).*(rm1)
        -27/16*(rp1).*(rp13).*(rm1)
        9/16*(rp1).*(rp13).*(rm13)
        ];
    end
    
    function ds = dShapeLine4( r )
      r=r(:)';
      rp1=r+1;
      rm1=r-1;
      rp13=r+1/3;
      rm13=r-1/3;
      ds = [
        -9/16*(rm13.*rp13 + rm1.*rp13 + rm1.*rm13)
        27/16*(rp1.*rm13 + rm1.*rp1 + rm1.*rm13)
        -27/16*(rp1.*rp13 + rm1.*rp1 + rm1.*rp13)
        9/16*(rp1.*rp13 + rp1.*rm13 + rm13.*rp13)
        ];
    end
    
  end % methods
  
  properties(Access=private)
    isOctave;
    intRule;
    numNodes, numElems, numElemNodes;
    numFEMEqns, totalNumEqns;
    numIntegrationPoints;
    hasODE, odeImpl;
    dirConsFlagsLeft, dirConsFlagsRight;
    y0FEM, y0;
    useDiagMassMat;
    vectorized;
    icDiagnostics, eqnDiagnostics, analyticalJacobian;
    pdeOpts;
    jacobianPattern, jacobianGroupList;
    sF, dSf; % shape function and derivative functions
    N, NN, dN; % shape functions and derivatives at int pts
    % temporary arrays for vectorized mode
    xPts;
  end
  
end % classdef

