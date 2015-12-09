## -*- texinfo -*-
## @deftypefn  {Function File} {@var{collist} = } monochromaticsoln (@var{col}, @var{eqn})
## @deftypefnx {Function File} {@var{collist} = } monochromaticsoln (@dots{}, 'Constant', @var{cons})
## @deftypefnx {Function File} {@var{collist} = } monochromaticsoln (@dots{}, 'Relator', @var{rel})
## @deftypefnx {Function File} {@var{collist} = } monochromaticsoln (@dots{}, 'verbose')
## Utilizes a backtracking algorithm to obtain a list of the valid
## monochromatic solution avoiding colorings of maximal length for
## the equation specified by @var{eqn}, @var{rel}, and @var{cons},
##
## @nospell{@var{eqn}(1)*@var{x}(1) + @var{eqn}(2)*@var{x}(2) + @dots{} + @var{eqn}(end)*@var{x}(end) + @var{cons} = 0}
##
## @end deftypefn
function solmat = monochromaticsoln ( varargin )
	p=inputParser();
	p.FunctionName='monochromaticsoln';
	p=p.addRequired('ColorNo',@isnumeric);
	p=p.addRequired('CoefVector',@isnumeric);
	p=p.addParamValue('Constant',0,@isnumeric);
	validate_Relator = @(x) any(strcmp(x,{'=','<','<='}));
	p=p.addParamValue('Relator','=',validate_Relator);
	p=p.addSwitch('verbose');
	p=p.parse(varargin{:});

	ColorNo=p.Results.ColorNo;
	CoefVector=p.Results.CoefVector;
	Constant=p.Results.Constant;
	Relator=p.Results.Relator;
	verbose=p.Results.verbose;

	% Validity check
	if( length(CoefVector) > 7)
		error('CoefVector too large!');
	endif
	% Relator is currently ignored.
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
		if( coloring(end)==0 )
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
