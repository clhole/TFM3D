function [ T ] = Traction( imgpair, Tparam )
%TRACTION     calculates the traction stresses from a given displacement field
%Input:
%  <imgpair>  image pair generated by <CreateBeadImagePair.m> including
%             displacement field data as calculated by <Displacement.m>
%  <Tparam>   structure containing the traction field calculation parameters:
%             .E               Young's modulus of the substrate [Pa]       [nx1]
%             .nu              Poisson's ratio                             [nx1]
%             .Talg            traction field calculation algorithms       {nx1}
%                              (options: 'direct','LsqFit','FEA')
%             .Direct.spacing  direct: grid spacing for grad(u) [px]       [nx1]
%             .LsqFit.winsize  LsqFit: window size for grad(u) [px]        [nx1]
%             .FEA.elemsize    FEA: surface element size [cell,gel][�m]    [nx2]
%             .Green.elemsize  Green: surface element size [cell,gel][�m]  [nx2]
%Output:
%  <T>        structure containing the traction stress field [pN/�m^2 = Pa]
%CL

% Get parameters
E = Tparam.E;
nu = Tparam.nu;
Talg = Tparam.Talg;
ualg = imgpair.uparam.ualg;
nualg = size(ualg,2);
imgparam = imgpair.imgparam;
range = imgparam.range;
fvcell = imgparam.fvcell;

% Determine mode for FEA stress calculation
if isfield(imgpair.imgparam,'fuFEA')
    sim = true;
else
    sim = false;
end

% Get labels and number of displacement evaluation methods and parameter sets
for ua = 1:nualg
    fieldlbl = imgpair.u.(ualg{ua}).Properties.VariableNames;
    ueval{ua} = fieldlbl(length(fieldnames(imgpair.uparam.(ualg{ua})))+1:end);
    nueval(ua) = size(ueval{ua},2);
    nuparamsets(ua) = size(imgpair.u.(ualg{ua}),1);
end

% Build the data structure using the displacement data structure as a template
T = imgpair.u;
% Delete displacement data
for ua = 1:nualg
    for ue = 1:nueval(ua)
        for up = 1:nuparamsets
            T.(ualg{ua}).(ueval{ua}{ue}){up} = [];
        end
    end
end

% Repeat for all displacement calculation algorithms
for ua = 1:nualg
    % Repeat for all displacement evaluation methods
    for ue = 1:nueval(ua)
        % Repeat for all displacement calculation parameter sets
        for up = 1:nuparamsets(ua)
            
            % Get positions and corresponding displacements
            pos = imgpair.u.(ualg{ua}).(ueval{ua}{ue}){up}.data{:}.pos;
            ucalc = imgpair.u.(ualg{ua}).(ueval{ua}{ue}){up}.data{:}.ucalc;
            
            % Remove data containing NaNs
            u = [pos,ucalc];
            u(any(isnan(u),2),:) = [];
            pos = u(:,1:3);
            ucalc = u(:,4:6);

            % Direct method 1: Direct differentiation of displacement field --> strain field --> traction field
            %****************************************************************************************************
            if any(strcmp('Direct',Talg))               
                
                T.(ualg{ua}).(ueval{ua}{ue}){up}.Direct = table();
                
                for i = 1:size(Tparam.Direct.spacing,2)
                    
                    spacing = Tparam.Direct.spacing(i);
                    
                    % Calculate traction field
                    fprintf('Traction field calculation: Direct with spacing = %.1f �m',spacing)
                    [pos_Direct,s_Direct] = DirectStressCalc(pos,ucalc,E,nu,spacing);

                    % Remove data points within the cell
                    inshp = inpolyhedron(fvcell,pos_Direct);
                    [row,~] = find(inshp);
                    pos_Direct(row,:) = [];
                    s_Direct(row,:) = [];
                    
                    % Evaluate traction stress field
                    mode = 'FDM';
                    evalDirect = TractionEval2(pos_Direct,s_Direct,imgparam,mode);

                    % Write data to table
                    T.(ualg{ua}).(ueval{ua}{ue}){up}.Direct = [T.(ualg{ua}).(ueval{ua}{ue}){up}.Direct;...
                        [table(spacing),evalDirect]];
                
                end
                
            end

            % Direct method 2: Local LSQ fitted displacement gradients --> strain field --> traction field 
            %****************************************************************************************************
            if any(strcmp('LsqFit',Talg))
                
                T.(ualg{ua}).(ueval{ua}{ue}){up}.LsqFit = table();
                
                for i = 1:size(Tparam.LsqFit.spacing,2)
                    for j = 1:size(Tparam.LsqFit.winsize,2)

                        % Get parameters
                        spacing = Tparam.LsqFit.spacing(i);
                        winsize = Tparam.LsqFit.winsize(j);

                        if strncmp('grid',ueval{ua}{ue},4) % LsqFit only works on gridded data

                            % Calculate traction field
                            fprintf('Traction field calculation: LsqFit with spacing = %.1f �m, winsize = %.1f\n',...
                                spacing, winsize)
                            [pos_LsqFit,s_LsqFit] = LsqFitStressCalc(pos,ucalc,E,nu,spacing,winsize);

                            % Remove data points within the cell
                            inshp = inpolyhedron(fvcell,pos_LsqFit);
                            [row,~] = find(inshp);
                            pos_LsqFit(row,:) = [];
                            s_LsqFit(row,:) = [];
                            
                            % Evaluate traction stress field
                            mode = 'LSQ';
                            evalLsqFit = TractionEval2(pos_LsqFit,s_LsqFit,imgparam,mode);

                            % Write data to table
                            T.(ualg{ua}).(ueval{ua}{ue}){up}.LsqFit = [T.(ualg{ua}).(ueval{ua}{ue}){up}.LsqFit;...
                                [table(spacing,winsize),evalLsqFit]];
                            
                        else

                            T.(ualg{ua}).(ueval{ua}{ue}){up}.LsqFit = [T.(ualg{ua}).(ueval{ua}{ue}){up}.LsqFit;...
                                [table(spacing,winsize),table({table(NaN)},{table(NaN)},'VariableNames',{'substrate','cell'})]];
                            
                        end

                    end
                end
                
            end                 
            
            % Direct method 3: FEA with displacement field as boundary condition --> traction field
            %****************************************************************************************************
            if any(strcmp('FEA',Talg))

                T.(ualg{ua}).(ueval{ua}{ue}){up}.FEA = table();

                for i = 1:size(Tparam.FEA.elemsize,1)
                    
                    % Get parameters
                    elemsize = Tparam.FEA.elemsize(i,:);

                    % Calculate traction field
                    fprintf('Traction field calculation: FEA with elemsize = [%.1f,%.1f]\n',...
                        elemsize(1),elemsize(2))
%                     f1 = scatteredInterpolant(pos,ucalc(:,1));
%                     f2 = scatteredInterpolant(pos,ucalc(:,2));
%                     f3 = scatteredInterpolant(pos,ucalc(:,3));
%                     ut(:,1) = f1(fvcell.vertices);
%                     ut(:,2) = f2(fvcell.vertices);
%                     ut(:,3) = f3(fvcell.vertices);
%                     fvcell.vertices = fvcell.vertices - ut;
                    [pos_FEA,s_FEA] = FEAStressCalc(pos,ucalc,E,nu,elemsize(1),elemsize(2),range,fvcell);

                    % Evaluate traction stress field
                    mode = 'FEA';
                    evalFEA = TractionEval2(pos_FEA,s_FEA,imgparam,mode);

                    % Write data to table
                    T.(ualg{ua}).(ueval{ua}{ue}){up}.FEA = [T.(ualg{ua}).(ueval{ua}{ue}){up}.FEA;...
                        [table(elemsize),evalFEA]];

                end

            end
            
            % Indirect method: Determine discretized Green's function and solve u=GT for cellular tractions
            %**************************************************************************************************
            if any(strcmp('Green',Talg))

                T.(ualg{ua}).(ueval{ua}{ue}){up}.Green = table();

                for i = 1:size(Tparam.Green.elemsize,1)
                    
                    % Get parameters
                    elemsize = Tparam.Green.elemsize(i,:);
                    

                    % Calculate traction field
                    fprintf('Traction field calculation: Green with elemsize = [%.1f,%.1f]\n',...
                        elemsize(1),elemsize(2))
                    [pos_Green,s_Green] = GreenStressCalc(pos,ucalc,E,nu,elemsize(1),elemsize(2),fvcell);

                    % Evaluate traction stress field
                    mode = 'GRE';
                    evalGreen = TractionEval2(pos_Green,s_Green,imgparam,mode);

                    % Write data to table
                    T.(ualg{ua}).(ueval{ua}{ue}){up}.Green = [T.(ualg{ua}).(ueval{ua}{ue}){up}.Green;...
                        [table(elemsize),evalGreen]];

                end

            end
            

        end %up
    end %ue
end %ua

end
