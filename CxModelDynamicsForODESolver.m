%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shaul Druckmann's code for the ODE of neural activity in the CX
% 05/03/2016: commented by DTE
% 05/13/2016: DTE added an input current term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = CxModelDynamicsForODESolver(t,y,fhandle,paramStruct,J,Iin)
  tau = paramStruct.tau;
  dy = -(y./tau)+(1./tau).*(J*fhandle(y,paramStruct))+Iin./tau;
end