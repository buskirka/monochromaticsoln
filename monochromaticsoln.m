## -*- texinfo -*-
## @deftypefn  {Function File} {@var{collist} = } monochromaticsoln (@var{col}, @var{eqn})
## @deftypefnx {Function File} {@var{collist} = } monochromaticsoln (@dots{}, 'Constant', @var{cons})
## @deftypefnx {Function File} {@var{collist} = } monochromaticsoln (@dots{}, 'Relator', @var{rel})
## @deftypefnx {Function File} {@var{collist} = } monochromaticsoln (@dots{}, 'verbose')
## Utilize a backtracking algorithm to obtain a list of the valid
## monochromatic solution avoiding colorings of maximal length for
## the equation specified by @var{eqn}, @var{rel}, and @var{cons},
##
## @nospell{@var{eqn}(1)*@var{x}(1) + @var{eqn}(2)*@var{x}(2) + @dots{} + @var{eqn}(end)*@var{x}(end) + @var{cons} = 0}
##
## @noindent
## where '=' may be substituted for the relator specified by @var{rel}.
##
## The @var{col} input may be one of two types of input:
## @itemize
## @item If @var{col} is an integer, it determines the number of colors we want in our
## colorings. The first @var{col} positive integers will then be our set of colors.
## @item If @var{col} is a @nospell{1xN} vector, it will be interpreted as a starting
## coloring for the procedure. In this case the color set will be taken to be
## @nospell{max(@var{col})}.
## @end itemize
##
## @var{eqn} is a @nospell{1xN} vector specifying the coefficients used in the equation 
## or inequality which we avoid monochromatic solutions to. These may be nonzero integer or
## floating point values, but using noninteger values may cause unpredictable behavior.
## @var{eqn} should also be relatively short, as otherwise the solution tables used
## will become excessively large.
##
## @var{cons} should be a numeric value, and most likely an integer.  
## Noninteger values are permitted but may cause unpredictable behavior.
##
## @var{rel} can be any of "=", "<", or "<=". To use the other inequalities, if the user
## desires, one may (of course) negate their @var{eqn} and @var{cons}. 
##
## Use of the 'verbose' flag causes the function to print to standard output a bit of 
## useful information for determining the function's progress towards its completion,
## and which can be used to partially recover from an abrupt cancellation. 
##
## @itemize
## @item To compute the 3rd Schur number, one may execute
## @example
## S(3)=monochromaticsoln(3,[1,1,-1])
## @end example
## @item To compute valid maximum-length MCS-avoiding 3-colorings to the inequality 
## @nospell{x1+x2+x3+1<x4} with
## verbose output,
## @example
## monochromaticsoln(3, [1,1,1,-1], ...
##                   'Constant', 1, ...
##                   'Relator', '<', ...
##                   'verbose')
## @end itemize
##
## @end deftypefn
function solmat = monochromaticsoln ( varargin )
	p=inputParser();
	p.FunctionName='monochromaticsoln';
	p.addRequired('ColorNo',@isnumeric);
	p.addRequired('CoefVector',@isnumeric);
	p.addParamValue('Constant',0,@isnumeric);
	validate_Relator = @(x) any(strcmp(x,{'=','<','<='}));
	p.addParamValue('Relator','=',validate_Relator);
	p.addSwitch('verbose');
	p.parse(varargin{:});

	ColorNo=p.Results.ColorNo;
	CoefVector=p.Results.CoefVector;
	Constant=p.Results.Constant;
	Relator=p.Results.Relator;
	verbose=p.Results.verbose;

	% Validity check
	if( length(CoefVector) > 7)
		error('CoefVector too large!');
	endif
	% The relator to be used.
    if(  (Relator~='=') && (Relator~='<') && (Relator~='<=') )
		warning('Erroneous Relator defaulted to =.');
		Relator='=';
	endif
	% Print the equation to the user.
	string=[num2str(CoefVector(1)),' x(1)'];
	for ind=2:length(CoefVector)
		val=CoefVector(ind);
		if(val>0)
			op='+';
		elseif(val<0)
			op='-';
			val=-val;
		else
			error('Zero coefficient!');
		endif
		string=[string,' ',op,' ',num2str(val),' x(',num2str(ind),')'];
	endfor
	if(verbose) 
		printf([string,' + ',num2str(Constant),' = 0 \n']);
	else
		warning ("off", "Octave:broadcast");	
	endif
	
	solnstruct=[];
	highscores=[];
	% Construct the initial navigation matrix
	if( length(size(ColorNo)) > 2 )
		error('ColorNo too large!');
	elseif( size(ColorNo,1)==1 && size(ColorNo,2) == 1 )
		navmat=( repmat((1:ColorNo)',1,ColorNo) <= repmat((1:ColorNo),ColorNo,1) );
	elseif( size(ColorNo,1)==1 )
		navmat=repmat((1:max(ColorNo))',1,length(ColorNo)) <= repmat(ColorNo,max(ColorNo),1);
		ColorNo=max(ColorNo);
	else
		navmat=(ColorNo==1);
		ColorNo=size(ColorNo,1);
	endif
	if(verbose)
		navmat
	endif
	
	tic
	counter=0;
	while(size(navmat,2)>1)
		% Construct the current coloring from the navmat.
		coloring=sum(navmat);
		% First we need to populate the solution matrix so that it at least
		% covers the coloring (and a bit more).
		if( max(size(solnstruct)) -1 < size(navmat,2) )
			clear solnstruct;
			l=size(navmat,2)+5;
			s=CoefVector(1)*(1:l)'; 
			for i=2:length(CoefVector); 
				s = s + CoefVector(i)*permute((1:l)',[i,2:(i-1),1,(i+1):length(CoefVector)]);
			endfor
			if(Relator=='=')
				solnstruct=(s+Constant==0);
			elseif(Relator=='<')
				solnstruct=(s+Constant<0);
			elseif(Relator=='<=')
				solnstruct=(s+Constant<=0);
			endif
		endif

		% Check whether the current coloring is terminal (has no further coloring).
		if( min(coloring)==0 )
			% If the last column is invalid, delete it.
			navmat=navmat(:,1:(end-1));
			delete=max(find(navmat(:,end)));
			navmat(delete,end)=0;
		else
			mcs=0;
			% Check whether the current coloring is valid
			for i=1:ColorNo;
				c=find(coloring==i);
				if(length(c)>0)
					% Construct the appropriate subarray of the solution structure.
					if( length(CoefVector) == 1 )
						colm=solnstruct(c);
					elseif( length(CoefVector) == 2)
						colm=solnstruct(c,c);
					elseif( length(CoefVector) == 3)
						colm=solnstruct(c,c,c);
					elseif( length(CoefVector) == 4)
						colm=solnstruct(c,c,c,c);
					elseif( length(CoefVector) == 5)
						colm=solnstruct(c,c,c,c,c);
					elseif( length(CoefVector) == 6)
						colm=solnstruct(c,c,c,c,c,c);
					elseif( length(CoefVector) == 7)
						colm=solnstruct(c,c,c,c,c,c,c);
					else
						error(['Invalid length',num2str(length(CoefVector))]);
					endif
					% Sum up the subarray; count the number of MCS found.
					while(length(colm)>1)
						colm=sum(colm);
					end
					mcs += colm;
				endif
			endfor
			if(mcs==0)	
				%It's valid! Hurray!
				if(length(coloring)>size(highscores,2))
					highscores=coloring;
					if(verbose)
						printf('\n');
						highscores
						toc
					endif
				elseif(length(coloring)==size(highscores,2))
					highscores=[highscores;coloring];
				end
				navmat=[navmat,(1:ColorNo)'<=1+max(coloring)];
			else
				% A solution must have been found. The most recent attempt 
				% was not a valid coloring.
				delete=max(find(navmat(:,end)));
				navmat(delete,end)=0;
			endif
		endif
		%printf([char(8)*ones(1,14),'(',num2str(sum(sum(navmat))),',',num2str(length(coloring)),'/',num2str(size(highscores,2)),')  ']);
		counter+=1;
		if(verbose && mod(counter,5000)==0)
			str=mat2str(coloring);
			printf([mat2str(size(highscores)),' : ',str,'\n']);
		endif
		fflush(stdout);
	endwhile
	solmat=highscores;
endfunction
