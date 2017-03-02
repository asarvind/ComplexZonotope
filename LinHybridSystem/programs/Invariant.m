classdef Invariant
    
    properties
        NumLocs
        map 
        input
        PtopeTemplate
        edges
        stay
        initial
        
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
            txt1 = sprintf('[%s, %s]*X_%s ',pritempnames{1},sectempnames{1},tag);
            txt2 = sprintf('== [%s, %s]',pritempnames{2},sectempnames{2});
            txt3 = sprintf('*diag([%s; (%s-%s)/2]);',scalename,upperbound,lowerbound);
            retval = [txt1,txt2,txt3];
        end
        
        function [ retval ] = TransferCenter( tag,pritemps,sectemps,centers,lowerbounds,upperbounds )
            txt1 = sprintf('[%s, %s]*y_%s == ',pritemps{1},sectemps{1},tag);
            txt2 = sprintf('%s-%s',centers{2},centers{1});
            txt3 = sprintf('+%s*(%s+%s)/2',sectemps{2},upperbounds{2},lowerbounds{2});
            txt4 = sprintf('-%s*(%s+%s)/2;',sectemps{1},upperbounds{1},lowerbounds{1});
            retval = [txt1,txt2,txt3,txt4];
        end
        
        function [ retval ] = ConstraintACZtope( tag,scalename,upperbound,lowerbound )
            txt1 = sprintf('sum(abs(X_%s),2)+abs(y_%s) <= ',tag,tag);
            txt2 = sprintf('[%s; (%s-%s)/2];',scalename,upperbound,lowerbound);
            retval = [txt1,txt2];
        end
        
    end
    
    methods
                
        function this = Invariant( this,varargin )
            this = this;
        end
        
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
            retval = sprintf( '\n%s\n%s\n%s',txt1,txt2,txt3 ); 
        end
        
        function [ this ] = CreateFile( this )
            delete(this.filename)
            fid = fopen( this.filename,'a' );
            
            % begin_cvx
            fprintf( fid,'cvx_begin\n' );
            
            
            
            fprintf( fid,'\ncvx_end\n');
            % end_cvx
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end