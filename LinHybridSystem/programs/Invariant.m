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
            txt1 = sprintf('[%s, %s]*X%s ',pritempnames{1},sectempnames{1},tag);
            txt2 = sprintf('== [%s, %s]',pritempnames{2},sectempnames{2});
            txt3 = sprintf('*diag([%s; (%s-%s)/2]);',scalename,upperbound,lowerbound);
            retval = [txt1,txt2,txt3];
        end
        
        function [ retval ] = TransferCenter( tag,pritemps,sectemps,centers,lowerbounds,upperbounds )
            txt1 = sprintf('[%s, %s]*y%s == ',pritemps{1},sectemps{1},tag);
            txt2 = sprintf('%s-%s',centers{2},centers{1});
            txt3 = sprintf('+%s*(%s+%s)/2',sectemps{2},upperbounds{2},lowerbounds{2});
            txt4 = sprintf('-%s*(%s+%s)/2;',sectemps{1},upperbounds{1},lowerbounds{1});
            retval = [txt1,txt2,txt3,txt4];
        end
        
        function [ retval ] = ConstraintACZtope( tag,scalename,upperbound,lowerbound )
            txt1 = sprintf('sum(abs(X%s),2)+abs(y%s) <= ',tag,tag);
            txt2 = sprintf('[%s; (%s-%s)/2];',scalename,upperbound,lowerbound);
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
        
        function [ this ] = Invariant( system,filename )
            this.system = system;
            this.filename = filename;
            
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
            fprintf( fid,'\n global b_locmaxx');
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
            fprintf( fid,'cvx_begin\n' );
            
            n = H.dim;
            
            % declare primary variables
            for i = 1:H.NumLocs
                clear {'m1','m2','m3','m4'}
                m1 = size(V{i},2);
                m2 = size(W{i},2);
                fprintf( fid,'variables s_%s(%s,1), c_%s(%s,1)',num2str(i),num2str(m1),num2str(i),num2str(n) ); 
                fprintf( fid,'variables l_%s(%s,1), u_%s(%s,1)',num2str(i),num2str(m2),num2str(i),num2str(m2) );
            end
            
            % declare auxillary variables for all locations
            for i = 1:H.NumLocs
                clear {'m1','m2','m3','m4'}
                m1 = size(V{i},2);
                m2 = size(W{i},2);
                m3 = size(H.input{i}.V,2);
                fprintf( fid,'\n variables X_loc%s(%s,%s), y_loc%s(%s,1)',num2str(i),num2str(m1+m2),num2str(m1+m2+m3),num2str(i),num2str(n) );
                fprintf( fid,'\n variables lpr_loc%s(%s,1), upr_loc%s(%s,1)',num2str(i),num2str(m2),num2str(i),num2str(m2) );
            end
            
            % declare auxillary variables for all edges
            for i = 1:numel(H.edges)
                edge = H.edges{i};
                clear {'m1','m2','m3','m4'}
                m1 = size(V{edge.loc1},2);
                m2 = size(W{edge.loc1},2);
                m3 = size(V{edge.loc2},2);
                m4 = size(W{edge.loc2},2);
                fprintf( fid,'\n variables X_edge%s(%s,%s), y_edge%s(%s,1)',num2str(i),num2str(m1+m2),num2str(m3+m4),num2str(i),num2str(n) );
                fprintf( fid,'\n variables lpr_edge%s(%s,1), upr_edge%s(%s,1)',num2str(i),num2str(m4),num2str(i),num2str(m4) );
            end
            
            % declare optimization variables
            fprintf( fid,'\n variables epsilon,lambda \n');
            
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
                acz2.PrimTempname = sprintf( '[H.A{%s}*V{%s}, H.input{%s}.V]',num2str(i),num2str(i),num2str(i) );
                acz1.SecTempname = sprintf( 'H.A{%s}*W{%s}',num2str(i),num2str(i) );
                acz1.centername = sprintf( 'H.A{%s}*c_%s + H.input{%s}.c',num2str(i),num2str(i),num2str(i) );
                acz1.scalename = sprintf( '[s_%s; H.input{%s}.s]',num2str(i),num2str(i) );
                acz1.lowerbound = sprintf( 'l_%s',num2str(i) );
                acz1.upperbound = sprintf( 'u_%s',num2str(i) );
                
                % stay constraint before transition
                fprintf( fid,'\n l_%s >= H.stay{%s}.l',num2str(i),num2str(i) );
                fprintf( fid,'\n u_%s <= H.stay{%s}.u',num2str(i),num2str(i) );
                
                % Overapproximation after transition
                tag = sprintf( '_loc%s',num2str(i) );
                txt = this.IncludeACZtope( tag,acz1,acz2 );
                fprintf( fid,'\n %s \n',txt);
                
                % Invariance after transition
                fprintf( fid,'(1-b_locmin{%s}).*upr_loc%s+v_locmin{%s}-u_%s <= epsilon',num2str(i),num2str(i),num2str(i),num2str(i) );
                fprintf( fid,'(1-b_locmax{%s}).*lpr_loc%s+v_locmax{%s}-l_%s >= -1*epsilon',num2str(i),num2str(i),num2str(i),num2str(i) );
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
                acz2.PrimTempname = sprintf( 'H.edge{%s}.reset*V{%s}',num2str(i),num2str(edge.loc1) );
                acz1.SecTempname = sprintf( 'H.edge{%s}.reset*W{%s}',num2str(i),num2str(edge.loc1) );
                acz1.centername = sprintf( 'H.edge{%s}.reset*c_%s',num2str(i),num2str(edge.loc1) );
                acz1.scalename = sprintf( '[s_%s]',num2str(edge.loc1) );
                acz1.lowerbound = sprintf( '(1-b_premax{%s}).*l_%s+v_premax{%s}',num2str(i),num2str(edge.loc1),num2str(i) );
                acz1.upperbound = sprintf( '(1-b_premin{%s}).*u_%s+v_premin{%s}',num2str(i),num2str(edge.loc1),num2str(i) );
                
                % Overapproximation after transition
                tag = sprintf( '_edge%s',num2str(i) );
                txt = this.IncludeACZtope( tag,acz1,acz2 );
                fprintf( fid,'\n %s \n',txt);
                
                % Invariance after transition
                fprintf( fid,'(1-b_postmin{%s}).*upr_edge%s+v_postmin{%s}-u_%s <= epsilon',num2str(i),num2str(i),num2str(i),num2str(edge.loc2) );
                fprintf( fid,'(1-b_postmax{%s}).*lpr_edge%s+v_postmax{%s}-l_%s >= -1*epsilon',num2str(i),num2str(i),num2str(i),num2str(edge.loc2) );
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
                acz1.PrimTempname = sprintf( 'H.initial{%s}.V',num2str(i) );
                acz2.SecTempname = sprintf( 'H.initial{%s}.W',num2str(i) );
                acz2.centername = sprintf( 'H.initial{%s}.c',num2str(i) );
                acz2.scalename = sprintf( 'H.initial{%s}.s',num2str(i) );
                acz2.lowerbound = sprintf( 'H.initial{%s}.l',num2str(i) );
                acz2.upperbound = sprintf( 'H.initial{%s}.u',num2str(i) );
                
                % Containment of initial set
                tag = sprintf( '_init%s',num2str(i) );
                txt = this.IncludeACZtope( tag,acz1,acz2 );
                fprintf( fid,'\n %s \n',txt);
            end
            
            % safety constraints
            for i = 1:H.NumLocs
                txt2 = sprintf( '[s_%s; (u_%s-l_%s)/2]',num2str(i),num2str(i),num2str(i) );
                txt4 = sprintf( '(u_%s+l_%s)/2',num2str(i),num2str(i) );
                txt1 = sprintf( 'abs(H.safety{%s}.T*[V{%s}, W{%s}]*%s)',num2str(i),num2str(i),num2str(i),txt2 );
                txt3 = sprintf( 'H.safety{%s}.T*(c_%s+%s)',num2str(i),num2str(i),txt4 );
                txt5 = sprintf( 'lambda*H.safety{%s}.d',num2str(i) );
                fprintf( fid,'%s+%s <= %s',txt1,txt3,txt5 );
            end
            
            fprintf( fid,'\ncvx_end\n');
            % end_cvx
            
            for i = 1:H.numLocs
                fprintf( fid,'scale{%s} = s_%s',num2str(i),num2str(i) );
                fprintf( fid,'center{%s} = c_%s',num2str(i),num2str(i) );
                fprintf( fid,'lowerbounds{%s} = l_%s',num2str(i),num2str(i) );
                fprintf( fid,'upperbounds{%s} = u_%s',num2str(i),num2str(i) );
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