classdef Invariant
    
    properties
        system
        %.NumLocs
        %.dim
        %.map 
        %.input
        %  .V
        %  .W
        %  .s
        %  .c
        %  .l
        %  .u
        %.Ptemplate
        %.edges
        %   .map 
        %   .input
        %.stay
        %  .l
        %  .u
        %.initial
        %  .V
        %  .W
        %  .s
        %  .c
        %  .l
        %  .u
        %.safe
        %  .T
        %  .d
        
        filename
        PrimTemplate
        SecTemplate
        center
        scale
        lbs
        ubs
        % .V
        % .W
        % .c
        % .s
        % .l
        % .u
    end
    
    methods (Static)
        
        function [ retval ] = TransferBounds( tag,pritempnames,sectempnames,scalename,lowerbound,upperbound )
            txt1 = sprintf('[(%s), (%s)]*((Xreal%s)+1i*(Xcomp%s)) ',pritempnames{1},sectempnames{1},tag,tag);
            txt2 = sprintf('== [(%s), (%s)]',pritempnames{2},sectempnames{2});
            txt3 = sprintf('*diag([(%s); ((%s)-(%s))/2]);',scalename,upperbound,lowerbound);
            retval = [txt1,txt2,txt3];
        end
        
        function [ retval ] = TransferCenter( tag,pritemps,sectemps,centers,lowerbounds,upperbounds )
            txt1 = sprintf('[(%s), (%s)]*((yreal%s)+1i*(ycomp%s)) == ',pritemps{1},sectemps{1},tag,tag);
            txt2 = sprintf('(%s)-(%s)',centers{2},centers{1});
            txt3 = sprintf('+(%s)*((%s)+(%s))/2',sectemps{2},upperbounds{2},lowerbounds{2});
            txt4 = sprintf('-(%s)*((%s)+(%s))/2;',sectemps{1},upperbounds{1},lowerbounds{1});
            retval = [txt1,txt2,txt3,txt4];
        end
        
        function [ retval ] = ConstraintACZtope( tag,scalename,upperbound,lowerbound )
            txt1 = sprintf('sum(abs((Xreal%s)+1i*(Xcomp%s)),2)+abs((yreal%s)+1i*(ycomp%s)) <= ',tag,tag,tag,tag);
            txt2 = sprintf('[(%s); ((%s)-(%s))/2];',scalename,upperbound,lowerbound);
            retval = [txt1,txt2];
        end
        
    end
    
    methods
        
        function this = AddConstraints( this,tag,acz1,acz2 )
            fid = fopen( this.filename,'a' );
            txt = this.IncludeACZtope( tag,acz1,acz2 );
            fprintf( fid,txt );
            fclose(fid);
        end
        
        function [ retval ] = IncludeACZtope( this,tag,acz1,acz2 )
            
            txt1 = this.TransferBounds( tag,{acz1.PrimTempname,acz2.PrimTempname},{acz1.SecTempname,acz2.SecTempname},acz2.scalename,...
                                   acz2.lowerbound,acz2.upperbound );
            txt2 = this.TransferCenter( tag,{acz1.PrimTempname,acz2.PrimTempname},{acz1.SecTempname,acz2.SecTempname},{acz1.centername,acz2.centername},...
                                   {acz1.lowerbound,acz2.lowerbound},{acz1.upperbound,acz2.upperbound} );
            txt3 = this.ConstraintACZtope( tag,acz1.scalename,acz1.upperbound,acz1.lowerbound );
            retval = sprintf( '%s\n%s\n%s',txt1,txt2,txt3 ); 
        end
        
        function this = SynTemplate(this,E)
            % Collect all eigenvectors
            if isempty(E)
                for i = 1:this.system.NumLocs
                    [new,~] = eig(this.system.map{i});
                    E = [E,new];
                end
                for i = 1:numel(this.system.edges)
                    [~,new] = eig(this.system.edges{i}.map);
                    E = [E,new];
                end
            end
            
            % Compute the templates in each location
            for i = 1:this.system.NumLocs
                this.SecTemplate{i} = pinv(this.system.Ptemplate{i});
                Orth = null(this.system.Ptemplate{i});
                this.PrimTemplate{i} = (Orth)*(Orth)'*E;
            end
        end
        
        function [ this ] = Invariant( system,E,filename,type )
            clear -global
            this.system = system;
            this.filename = filename;
            
            % Synthesize template
            this = this.SynTemplate(E);
            
            % Create CVX file
            delete(this.filename)
            fid = fopen( this.filename,'a' );  
            
            global H;
            fprintf( fid,'\n global H');
            global V;
            fprintf( fid,'\n global V');
            global W;
            fprintf( fid,'\n global W');
            H = this.system;
            V = this.PrimTemplate;
            W = this.SecTemplate;
            
            global scale;
            fprintf( fid,'\n global scale');
            global center;
            fprintf( fid,'\n global center');
            global lowerbounds;
            fprintf( fid,'\n global lowerbounds');
            global upperbounds;
            fprintf( fid,'\n global upperbounds');
            
            % compute approxbounds for location dynamics
            global b_locmin;
            fprintf( fid,'\n global b_locmin');
            global b_locmax;
            fprintf( fid,'\n global b_locmax');
            global v_locmin;
            fprintf( fid,'\n global v_locmin');
            global v_locmax;
            fprintf( fid,'\n global v_locmax');
            for i = 1:H.NumLocs                
                [val1,val2] = init_minvect(H.stay{i}.u);
                [val3,val4] = init_maxvect(H.stay{i}.l);
                b_locmin{i} = val1; v_locmin{i} = val2;
                b_locmax{i} = val3; v_locmax{i} = val4;
            end
            
            
            % compute approxbounds for all edge transitions
            global bedge_premin;
            fprintf( fid,'\n global bedge_premin');
            global bedge_premax;
            fprintf( fid,'\n global bedge_premax');
            global bedge_postmin;
            fprintf( fid,'\n global bedge_postmin');
            global bedge_postmax;
            fprintf( fid,'\n global bedge_postmax');
            global vedge_premin;
            fprintf( fid,'\n global vedge_premin');
            global vedge_premax;
            fprintf( fid,'\n global vedge_premax');
            global vedge_postmin;
            fprintf( fid,'\n global vedge_postmin');
            global vedge_postmax;
            fprintf( fid,'\n global vedge_postmax');
            for i = 1:numel(H.edges)
                edge = H.edges{i};
                [val1,val2] = init_minvect(min(edge.u,H.stay{edge.loc1}.u));
                [val3,val4] = init_maxvect(max(edge.l,H.stay{edge.loc1}.l));
                [val5,val6] = init_minvect(H.stay{edge.loc2}.u);
                [val7,val8] = init_maxvect(H.stay{edge.loc2}.l);
                bedge_premin{i} = val1; vedge_premin{i} = val2;
                bedge_premax{i} = val3; vedge_premax{i} = val4;
                bedge_postmin{i} = val5; vedge_postmin{i} = val6;
                bedge_postmax{i} = val7; vedge_postmax{i} = val8;
            end
            
                      
     
            % begin_cvx
            fprintf( fid,'\n \n cvx_begin\n' );
            
            n = H.dim;
            
            comment = '%';
            
            % declare primary variables
            for i = 1:H.NumLocs
                clear {'m1','m2','m3','m4'}
                m1 = size(V{i},2);
                m2 = size(W{i},2);
                fprintf(fid,'\n %s declare primary variables for location %s',comment,num2str(i));
                fprintf( fid,'\n variables s_%s(%s,1) c_%s(%s,1)',num2str(i),num2str(m1),num2str(i),num2str(n) ); 
                fprintf( fid,'\n variables l_%s(%s,1) u_%s(%s,1)',num2str(i),num2str(m2),num2str(i),num2str(m2) );
                fprintf(fid,'\n');
            end
            
            % declare auxillary variables for all locations
            for i = 1:H.NumLocs
                clear {'m1','m2','m3','m4','m5','m6'}
                m1 = size(V{i},2);
                m2 = size(W{i},2);
                m3 = size(H.input{i}.V,2);
                if ~isempty(system.initial{i})
                    m4 = size(H.initial{i}.V,2);
                    m5 = size(H.initial{i}.W,2);
                end
                m6 = size(H.input{i}.W,2);
                fprintf(fid,'\n %s declare auxillary variables for location %s',comment,num2str(i));
                fprintf( fid,'\n variables Xreal_loc%s(%s,%s) yreal_loc%s(%s,1)',num2str(i),num2str(m1+m2),num2str(m1+m2+m3+m6),num2str(i),num2str(m1+m2) );
                fprintf( fid,'\n variables Xcomp_loc%s(%s,%s) ycomp_loc%s(%s,1)',num2str(i),num2str(m1+m2),num2str(m1+m2+m3+m6),num2str(i),num2str(m1+m2) );
                fprintf( fid,'\n variables Xreal_init%s(%s,%s) yreal_init%s(%s,1)',num2str(i),num2str(m1+m2),num2str(m4+m5),num2str(i),num2str(m1+m2) );
                fprintf( fid,'\n variables Xcomp_init%s(%s,%s) ycomp_init%s(%s,1)',num2str(i),num2str(m1+m2),num2str(m4+m5),num2str(i),num2str(m1+m2) );
                fprintf( fid,'\n variables lpr_loc%s(%s,1) upr_loc%s(%s,1)',num2str(i),num2str(m2),num2str(i),num2str(m2) );
                fprintf(fid,'\n');
            end
            
            % declare auxillary variables for all edges
            for i = 1:numel(H.edges)
                edge = H.edges{i};
                clear {'m1','m2','m3','m4'}
                m1 = size(V{edge.loc1},2);
                m2 = size(W{edge.loc1},2);
                m3 = size(V{edge.loc2},2);
                m4 = size(W{edge.loc2},2);
                m5 = size(H.edges{i}.input.V,2);
                m6 = size(H.edges{i}.input.W,2);
                fprintf(fid,'\n %s declare auxillary variables for edge %s',comment,num2str(i));
                fprintf( fid,'\n variables Xreal_edge%s(%s,%s) yreal_edge%s(%s,1)',num2str(i),num2str(m3+m4),num2str(m1+m2+m5+m6),num2str(i),num2str(m1+m2) );
                fprintf( fid,'\n variables Xcomp_edge%s(%s,%s) ycomp_edge%s(%s,1)',num2str(i),num2str(m3+m4),num2str(m1+m2+m5+m6),num2str(i),num2str(m1+m2) );
                fprintf( fid,'\n variables lpr_edge%s(%s,1) upr_edge%s(%s,1)',num2str(i),num2str(m4),num2str(i),num2str(m4) );
                fprintf(fid,'\n');
            end
            
            % declare optimization variables
            fprintf( fid,'\n variables epsilon lambda \n');
            
            % declare objective function
            switch type
                case 'optimize'
                    fprintf( fid,'\n minimize(lambda)\n     subject to\n');
                case 'feasibility'
                    fprintf( fid,'\n lambda == 1' );
            end
            
            % epsilon is not greater than zero
            fprintf( fid,'\n epsilon <= 0');
            
            % specify inclusion relations for intralocation dynamics
            for i = 1:H.NumLocs
                clear {'acz1','acz2','tag','txt'};
                
                % specify acz1 names
                acz1.PrimTempname = sprintf( 'V{%s}',num2str(i) );
                acz1.SecTempname = sprintf( 'W{%s}',num2str(i) );
                acz1.centername = sprintf( 'c_%s',num2str(i) );
                acz1.scalename = sprintf( 's_%s',num2str(i) );
                acz1.lowerbound = sprintf( 'lpr_loc%s',num2str(i) );
                acz1.upperbound = sprintf( 'upr_loc%s',num2str(i) );
                
                % specify acz2 names
                acz2.PrimTempname = sprintf( '[H.map{%s}*V{%s}, H.input{%s}.V]',num2str(i),num2str(i),num2str(i) );
                acz2.SecTempname = sprintf( '[H.map{%s}*W{%s}, H.input{%s}.W]',num2str(i),num2str(i),num2str(i) );
                acz2.centername = sprintf( 'H.map{%s}*c_%s + H.input{%s}.c',num2str(i),num2str(i),num2str(i) );
                acz2.scalename = sprintf( '[s_%s; H.input{%s}.s]',num2str(i),num2str(i) );
                acz2.lowerbound = sprintf( '[(1-b_locmax{%s}).*l_%s+v_locmax{%s}; H.input{%s}.l]',num2str(i),num2str(i),num2str(i),num2str(i) );
                acz2.upperbound = sprintf( '[(1-b_locmin{%s}).*u_%s+v_locmin{%s}; H.input{%s}.u]',num2str(i),num2str(i),num2str(i),num2str(i) );
                
                % Overapproximation after transition
                fprintf(fid,'\n %s overapproximation of location %s reach set',comment,num2str(i));
                tag = sprintf( '_loc%s',num2str(i) );
                txt = this.IncludeACZtope( tag,acz1,acz2 );
                fprintf( fid,'\n %s \n',txt);
                fprintf(fid,'\n');
                
                % Invariance after transition
                fprintf(fid,'\n %s invariance after transition for location %s',comment, num2str(i));
                fprintf( fid,'\n (1-b_locmin{%s}).*upr_loc%s+v_locmin{%s}-epsilon <= u_%s',num2str(i),num2str(i),num2str(i),num2str(i) );
                fprintf( fid,'\n l_%s-epsilon <= (1-b_locmax{%s}).*lpr_loc%s+v_locmax{%s}',num2str(i),num2str(i),num2str(i),num2str(i) );
                fprintf(fid,'\n');
            end
            
            % specify inclusion relations for interlocation dynamics
            for i = 1:numel(H.edges)
                edge = H.edges{i};
                
                clear {'acz1','acz2','tag','txt'};
                % specify acz1 names
                acz1.PrimTempname = sprintf( 'V{%s}',num2str(edge.loc2) );
                acz1.SecTempname = sprintf( 'W{%s}',num2str(edge.loc2) );
                acz1.centername = sprintf( 'c_%s',num2str(edge.loc2) );
                acz1.scalename = sprintf( 's_%s',num2str(edge.loc2) );
                acz1.lowerbound = sprintf( 'lpr_edge%s',num2str(i) );
                acz1.upperbound = sprintf( 'upr_edge%s',num2str(i) );
                
                % specify acz2 names
                acz2.PrimTempname = sprintf( '[H.edges{%s}.map*V{%s}, H.edges{%s}.input.V]',num2str(i),num2str(edge.loc1),num2str(i) );
                acz2.SecTempname = sprintf( '[H.edges{%s}.map*W{%s}, H.edges{%s}.input.W]',num2str(i),num2str(edge.loc1),num2str(i) );
                acz2.centername = sprintf( '(H.edges{%s}.map*c_%s+H.edges{%s}.input.c)',num2str(i),num2str(edge.loc1),num2str(i) );
                acz2.scalename = sprintf( '[s_%s; H.edges{%s}.input.s]',num2str(edge.loc1),num2str(i) );
                acz2.lowerbound = sprintf( '[(1-bedge_premax{%s}).*l_%s+vedge_premax{%s}; H.edges{%s}.input.l]',num2str(i),num2str(edge.loc1),num2str(i),num2str(i) );
                acz2.upperbound = sprintf( '[(1-bedge_premin{%s}).*u_%s+vedge_premin{%s}; H.edges{%s}.input.u]',num2str(i),num2str(edge.loc1),num2str(i),num2str(i) );
                
                % Overapproximation after transition
                fprintf(fid,'\n %s Overapproximation after edge %s transition',comment,num2str(i));
                tag = sprintf( '_edge%s',num2str(i) );
                txt = this.IncludeACZtope( tag,acz1,acz2 );
                fprintf( fid,'\n %s \n',txt);
                
                % Invariance after transition
                fprintf(fid,'\n %s Invariance condition after edge %s transition',comment,num2str(i));
                fprintf( fid,'\n (1-bedge_postmin{%s}).*upr_edge%s+vedge_postmin{%s}-epsilon <= u_%s ',num2str(i),num2str(i),num2str(i),num2str(edge.loc2) );
                fprintf( fid,'\n l_%s-epsilon <= (1-bedge_postmax{%s}).*lpr_edge%s+vedge_postmax{%s}',num2str(edge.loc2),num2str(i),num2str(i),num2str(i) );
                fprintf(fid,'\n');
            end
            
            % inclusion of initial set
            for i = 1:H.NumLocs
                clear {'acz1','acz2','tag','txt'};
                if ~isempty(system.initial{i})
                    % specify acz1 names
                    acz1.PrimTempname = sprintf( 'V{%s}',num2str(i) );
                    acz1.SecTempname = sprintf( 'W{%s}',num2str(i) );
                    acz1.centername = sprintf( 'c_%s',num2str(i) );
                    acz1.scalename = sprintf( 's_%s',num2str(i) );
                    acz1.lowerbound = sprintf( 'l_%s',num2str(i) );
                    acz1.upperbound = sprintf( 'u_%s',num2str(i) );

                    %specify acz2 names
                    acz2.PrimTempname = sprintf( 'H.initial{%s}.V',num2str(i) );
                    acz2.SecTempname = sprintf( 'H.initial{%s}.W',num2str(i) );
                    acz2.centername = sprintf( 'H.initial{%s}.c',num2str(i) );
                    acz2.scalename = sprintf( 'H.initial{%s}.s',num2str(i) );
                    acz2.lowerbound = sprintf( 'H.initial{%s}.l',num2str(i) );
                    acz2.upperbound = sprintf( 'H.initial{%s}.u',num2str(i) );

                    % Containment of initial set
                    fprintf(fid,'\n %s Containment of initial set in location %s',comment,num2str(i));
                    tag = sprintf( '_init%s',num2str(i) );
                    txt = this.IncludeACZtope( tag,acz1,acz2 );
                    fprintf( fid,'\n %s \n',txt);
                    fprintf(fid,'\n');
                end
            end
            
            % safety constraints
            for i = 1:H.NumLocs
                fprintf(fid,'\n %s safety condition in location %s',comment,num2str(i));
                txt2 = sprintf( '[s_%s; (u_%s-l_%s)/2]',num2str(i),num2str(i),num2str(i) );
                txt4 = sprintf( '(W{%s}*(u_%s+l_%s)/2)',num2str(i),num2str(i),num2str(i) );
                txt1 = sprintf( 'abs(H.safe{%s}.T*[V{%s}, W{%s}]*%s)',num2str(i),num2str(i),num2str(i),txt2 );
                txt3 = sprintf( 'H.safe{%s}.T*(c_%s+%s)',num2str(i),num2str(i),txt4 );
                txt5 = sprintf( 'lambda*H.safe{%s}.d',num2str(i) );
                fprintf( fid,'\n %s+%s <= %s\n',txt1,txt3,txt5 );
                fprintf(fid,'\n');
            end
            
            fprintf( fid,'\n cvx_end \n');
            % end_cvx
            
            for i = 1:H.NumLocs
                fprintf(fid,'\n %s assign optimized variables to ACZ in location %s', comment,num2str(i));
                fprintf( fid,'\n scale{%s} = s_%s;',num2str(i),num2str(i) );
                fprintf( fid,'\n center{%s} = c_%s;',num2str(i),num2str(i) );
                fprintf( fid,'\n lowerbounds{%s} = l_%s;',num2str(i),num2str(i) );
                fprintf( fid,'\n upperbounds{%s} = u_%s;',num2str(i),num2str(i) );
                fprintf(fid,'\n');
            end            
            fclose(fid);
            
            run(this.filename);
                     
            this.scale = scale;
            this.center = center;
            this.lbs = lowerbounds;
            this.ubs = upperbounds;
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end