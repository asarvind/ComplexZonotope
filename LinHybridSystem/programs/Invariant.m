classdef Invariant
    
    properties
        system
        %.NumLocs
        %.dim
        %.map 
        %.input
        %.Ptemplate
        %.edges
        %.stay
        %.initial
        %.safety
        
        filename
        
        PrimTemplate
        SecTemplate
        center
        scale
        lbs
        ubs
    end
    
    methods (Static)
        
        function [ retval ] = TransferBounds( tag,pritempnames,sectempnames,scalename,lowerbound,upperbound )
            txt1 = sprintf('[(%s), (%s)]*(Xreal%s+1i*Xcomp%s) ',pritempnames{1},sectempnames{1},tag,tag);
            txt2 = sprintf('== [(%s), (%s)]',pritempnames{2},sectempnames{2});
            txt3 = sprintf('*diag([(%s); ((%s)-(%s))/2]);',scalename,upperbound,lowerbound);
            retval = [txt1,txt2,txt3];
        end
        
        function [ retval ] = TransferCenter( tag,pritemps,sectemps,centers,lowerbounds,upperbounds )
            txt1 = sprintf('[(%s), (%s)]*(yreal%s+1i*ycomp%s) == ',pritemps{1},sectemps{1},tag,tag);
            txt2 = sprintf('(%s)-(%s)',centers{2},centers{1});
            txt3 = sprintf('+(%s)*((%s)+(%s))/2',sectemps{2},upperbounds{2},lowerbounds{2});
            txt4 = sprintf('-(%s)*((%s)+(%s))/2;',sectemps{1},upperbounds{1},lowerbounds{1});
            retval = [txt1,txt2,txt3,txt4];
        end
        
        function [ retval ] = ConstraintACZtope( tag,scalename,upperbound,lowerbound )
            txt1 = sprintf('sum(abs(Xreal%s+1i*Xcomp%s),2)+abs(yreal%s+1i*ycomp%s) <= ',tag,tag,tag,tag);
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
        
        function this = SynTemplate(this)
            % Collect all eigenvectors
            E = [];
            for i = 1:this.system.NumLocs
                [~,new] = eig(this.system.map{i});
                E = [E,new];
            end
            for i = 1:numel(this.system.edges)
                [~,new] = eig(this.system.edges{i}.reset);
                E = [E,new];
            end
            
            % Compute the templates in each location
            for i = 1:this.system.NumLocs
                this.SecTemplate{i} = pinv(this.system.Ptemplate{i});
                Orth = null(this.system.Ptemplate{i});
                this.PrimTemplate{i} = (Orth)*(Orth)'*E;
            end
        end
        
        function [ this ] = Invariant( system,filename )
            this.system = system;
            this.filename = filename;
            
            % Synthesize template
            this = this.SynTemplate();
            
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
                [b_locmin,v_locmin] = init_minvect(H.stay{i}.u);
                [b_locmax,v_locmax] = init_maxvect(H.stay{i}.l);
            end
            
            
            % compute approxbounds for all edge transitions
            global b_premin;
            fprintf( fid,'\n global b_premin');
            global b_premax;
            fprintf( fid,'\n global b_premax');
            global b_postmin;
            fprintf( fid,'\n global b_postmin');
            global b_postmax;
            fprintf( fid,'\n global b_postmax');
            global v_premin;
            fprintf( fid,'\n global v_premin');
            global v_premax;
            fprintf( fid,'\n global v_premax');
            global v_postmin;
            fprintf( fid,'\n global v_postmin');
            global v_postmax;
            fprintf( fid,'\n global v_postmax');
            for i = 1:numel(H.edges)
                edge = H.edges{i};
                [b_premin{i},v_premin{i}] = init_minvect(min(edge.u,H.stay{edge.loc1}.u));
                [b_premax{i},v_premax{i}] = init_maxvect(max(edge.l,H.stay{edge.loc1}.l));
                [b_postmin{i},v_postmin{i}] = init_minvect(H.stay{edge.loc2}.u);
                [b_postmax{i},v_postmax{i}] = init_maxvect(H.stay{edge.loc2}.l);
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
                clear {'m1','m2','m3','m4'}
                m1 = size(V{i},2);
                m2 = size(W{i},2);
                m3 = size(H.input{i}.V,2);
                fprintf(fid,'\n %s declare auxillary variables for location %s',comment,num2str(i));
                fprintf( fid,'\n variables Xreal_loc%s(%s,%s) yreal_loc%s(%s,1)',num2str(i),num2str(m1+m2),num2str(m1+m2+m3),num2str(i),num2str(n) );
                fprintf( fid,'\n variables Xcomp_loc%s(%s,%s) ycomp_loc%s(%s,1)',num2str(i),num2str(m1+m2),num2str(m1+m2+m3),num2str(i),num2str(n) );
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
                fprintf(fid,'\n %s declare auxillary variables for edge %s',comment,num2str(i));
                fprintf( fid,'\n variables Xreal_edge%s(%s,%s) yreal_edge%s(%s,1)',num2str(i),num2str(m1+m2),num2str(m3+m4),num2str(i),num2str(n) );
                fprintf( fid,'\n variables Xcomp_edge%s(%s,%s) ycomp_edge%s(%s,1)',num2str(i),num2str(m1+m2),num2str(m3+m4),num2str(i),num2str(n) );
                fprintf( fid,'\n variables lpr_edge%s(%s,1) upr_edge%s(%s,1)',num2str(i),num2str(m4),num2str(i),num2str(m4) );
                fprintf(fid,'\n');
            end
            
            % declare optimization variables
            fprintf( fid,'\n variables epsilon lambda \n');
            
            % declare objective function
            fprintf( fid,'\n minimize(lambda)\n     subject to\n');
            
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
                acz2.SecTempname = sprintf( 'H.map{%s}*W{%s}',num2str(i),num2str(i) );
                acz2.centername = sprintf( 'H.map{%s}*c_%s + H.input{%s}.c',num2str(i),num2str(i),num2str(i) );
                acz2.scalename = sprintf( '[s_%s; H.input{%s}.s]',num2str(i),num2str(i) );
                acz2.lowerbound = sprintf( 'l_%s',num2str(i) );
                acz2.upperbound = sprintf( 'u_%s',num2str(i) );
                
                % stay constraint before transition
                fprintf(fid,'\n %s stay constraint before transition for location %s',comment,num2str(i));
                fprintf( fid,'\n l_%s >= H.stay{%s}.l',num2str(i),num2str(i) );
                fprintf( fid,'\n u_%s <= H.stay{%s}.u',num2str(i),num2str(i) );
                fprintf(fid,'\n');
                
                % Overapproximation after transition
                fprintf(fid,'\n %s overapproximation of location %s reach set',comment,num2str(i));
                tag = sprintf( '_loc%s',num2str(i) );
                txt = this.IncludeACZtope( tag,acz1,acz2 );
                fprintf( fid,'\n %s \n',txt);
                fprintf(fid,'\n');
                
                % Invariance after transition
                fprintf(fid,'\n %s invariance after transition for location %s',comment, num2str(i));
                fprintf( fid,'\n (1-b_locmin{%s}).*upr_loc%s+v_locmin{%s}-u_%s <= epsilon',num2str(i),num2str(i),num2str(i),num2str(i) );
                fprintf( fid,'\n (1-b_locmax{%s}).*lpr_loc%s+v_locmax{%s}-l_%s >= -1*epsilon',num2str(i),num2str(i),num2str(i),num2str(i) );
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
                acz2.PrimTempname = sprintf( 'H.edges{%s}.reset*V{%s}',num2str(i),num2str(edge.loc1) );
                acz2.SecTempname = sprintf( 'H.edges{%s}.reset*W{%s}',num2str(i),num2str(edge.loc1) );
                acz2.centername = sprintf( 'H.edges{%s}.reset*c_%s',num2str(i),num2str(edge.loc1) );
                acz2.scalename = sprintf( 's_%s',num2str(edge.loc1) );
                acz2.lowerbound = sprintf( '(1-b_premax{%s}).*l_%s+v_premax{%s}',num2str(i),num2str(edge.loc1),num2str(i) );
                acz2.upperbound = sprintf( '(1-b_premin{%s}).*u_%s+v_premin{%s}',num2str(i),num2str(edge.loc1),num2str(i) );
                
                % Overapproximation after transition
                fprintf(fid,'\n %s Overapproximation after edge %s transition',comment,num2str(i));
                tag = sprintf( '_edge%s',num2str(i) );
                txt = this.IncludeACZtope( tag,acz1,acz2 );
                fprintf( fid,'\n %s \n',txt);
                
                % Invariance after transition
                fprintf(fid,'\n %s Invariance condition after edge %s transition',comment,num2str(i));
                fprintf( fid,'\n (1-b_postmin{%s}).*upr_edge%s+v_postmin{%s}-u_%s <= epsilon',num2str(i),num2str(i),num2str(i),num2str(edge.loc2) );
                fprintf( fid,'\n (1-b_postmax{%s}).*lpr_edge%s+v_postmax{%s}-l_%s >= -1*epsilon',num2str(i),num2str(i),num2str(i),num2str(edge.loc2) );
                fprintf(fid,'\n');
            end
            
            % inclusion of initial set
            for i = 1:H.NumLocs
                clear {'acz1','acz2','tag','txt'};
                
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
            
            % safety constraints
            for i = 1:H.NumLocs
                fprintf(fid,'\n %s safety condition in location %s',comment,num2str(i));
                txt2 = sprintf( '[s_%s; (u_%s-l_%s)/2]',num2str(i),num2str(i),num2str(i) );
                txt4 = sprintf( '(u_%s+l_%s)/2',num2str(i),num2str(i) );
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
            
%             run(this.filename);
%                      
%             this.scale = scale;
%             this.center = center;
%             this.lbs = lowerbounds;
%             this.ubs = upperbounds;
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end