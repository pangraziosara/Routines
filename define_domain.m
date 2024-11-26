

function domain_structure=define_domain(domain_str,varargin)

domain_structure.domain=domain_str;

args = varargin;
L=length(args);
k=1;

switch domain_str

    case 'polygon'

        % ........................ setting defaults .......................
        %
        % domain_structure.vertices:
        % 0. it is a matrix N x 2;
        % 1. each row consists of the coordinates of a vertex of the domain
        %    in cartesian coordinates;

        % setting defaults
        example_type=6;
        domain_structure=gallery_domains_2D(domain_str,example_type);

        while k <= L
            strL=args{k};
            switch strL
                case 'vertices'
                    domain_structure.vertices=args{k+1};
            end
            k=k+2;
        end

        vertices=domain_structure.vertices;

        if norm(vertices(1,:)-vertices(end,:)) == 0
            vertices=vertices(1:end-1,:);
        end

        XV=vertices(:,1); YV=vertices(:,2);

        domain_structure.polyshape=polyshape(XV,YV);



    case 'disk'

        % ........................ setting defaults .......................
        %
        % domain_structure.center:
        % 0. it is a vector 1 x 2; it consists of the coordinates of the
        % center of the circular annular sector;
        %
        % domain_structure.radius:
        % 0. it is a vector of length 1;
        % 1. first value consists of the radius of the first disk;
        % 2. second value consists of the radius of the second disk;


        % setting defaults
        domain_structure.center=[0 0];
        domain_structure.radius=1;

        while k <= L
            strL=args{k};
            switch strL
                case 'center'
                    % setting center as a row vector
                    value=args{k+1}; value=(value(:))';
                    domain_structure.center=value;
                case 'radius'
                    domain_structure.radius=args{k+1};
            end
            k=k+2;
        end


        % ellblend structure (as sector)

        R2=0;
        R1=domain_structure.radius;
        C=domain_structure.center;

        domain_structure.A=[R1 0; R2 0];
        domain_structure.B=[0 R1; 0 R2];
        domain_structure.C=[C; C];
        domain_structure.alpha=0;
        domain_structure.beta=2*pi;



    case 'lune'

        % ........................ setting defaults .......................
        %
        % domain_structure.centers:
        % 0. it is a matrix 2 x 2;
        % 1. first row consists of the coordinates of the center of the
        %    first lune;
        % 2. second row consists of the coordinates of the center of the
        %    second lune;
        %
        % domain_structure.radii:
        % 0. it is a vector of length 2;
        % 1. first value consists of the radius of the first lune;
        % 2. second value consists of the radius of the second lune;

        domain_structure.centers=[0.5 0.5; 1.1 0.5];
        domain_structure.radii=[0.5; 0.5];

        while k <= L
            strL=args{k};
            switch strL
                case 'centers'
                    domain_structure.center=args{k+1};
                case 'radii'
                    domain_structure.radius=args{k+1};
            end
            k=k+2;
        end

        
        domain_structure.NURBS=ellblend2NURBS(domain_structure);
        


    case 'circular-annular-sector'

        % ........................ setting defaults .......................
        %
        % domain_structure.center:
        % 0. it is a vector 1 x 2; it consists of the coordinates of the
        % center of the circular annular sector;
        %
        % domain_structure.radii:
        % 0. it is a vector of length 2;
        % 1. first value consists of the radius of the first disk;
        % 2. second value consists of the radius of the second disk;
        %
        % domain_structure.angles:
        % 0. it is a vector of length 2;
        % 1. the first component is the parameter "alpha" and must be in
        %    [-pi,pi]
        % 2. the second component is the parameter "beta" and must be in
        %    [-pi,pi]


        % default example:
        domain_structure.center=[0 0];
        domain_structure.radii=[0.5 2];
        domain_structure.angles=[pi/4; pi/2];

        while k <= L
            strL=args{k};
            switch strL
                case 'center'
                    % setting center as a row vector
                    value=args{k+1}; value=(value(:))';
                    domain_structure.center=value;
                case 'radii'
                    domain_structure.radius=args{k+1};
                case 'angles'
                    domain_structure.angles=args{k+1};
            end
            k=k+2;
        end

        domain_structure.NURBS=ellblend2NURBS(domain_structure);


        % ellblend structure

        R=domain_structure.radii;
        R1=R(1);
        R2=R(2);
        C=domain_structure.center;
        thetaV=domain_structure.angles;

        domain_structure.A=[R1 0; R2 0];
        domain_structure.B=[0 R1; 0 R2];
        domain_structure.C=[C; C];
        domain_structure.alpha=thetaV(1);
        domain_structure.beta=thetaV(2);




    case 'sector'

        % ........................ setting defaults .......................
        %
        % domain_structure.center:
        % 0. it is a vector 1 x 2; it consists of the coordinates of the
        % center of the sector;
        %
        % domain_structure.radius:
        % 0. it is a vector of length 1;
        % 1. it consists of the radius of the disk;
        %
        % domain_structure.angles:
        % 0. it is a vector of length 2;
        % 1. the first component is the parameter "alpha" and must be in
        %    [-pi,pi]
        % 2. the second component is the parameter "beta" and must be in
        %    [-pi,pi]

        % default example:
        domain_structure.center=[1 1.5];
        domain_structure.radius=1;
        domain_structure.angles=[pi/4; pi/2];

        while k <= L
            strL=args{k};
            switch strL
                case 'center'
                    % setting center as a row vector
                    value=args{k+1}; value=(value(:))';
                    domain_structure.center=value;
                case 'radius'
                    domain_structure.radius=args{k+1};
                case 'angles'
                    domain_structure.angles=args{k+1};
            end
            k=k+2;
        end

        domain_structure.NURBS=ellblend2NURBS(domain_structure);


        % ellblend structure (as symmetric annular sector)

        R1=0;
        R2=domain_structure.radius;
        C=domain_structure.center;
        thetaV=domain_structure.angles;

        domain_structure.A=[R1 0; R2 0];
        domain_structure.B=[0 R1; 0 R2];
        domain_structure.C=[C; C];
        domain_structure.alpha=thetaV(1);
        domain_structure.beta=thetaV(2);




    case 'asymmetric-circular-sector'

        % ........................ setting defaults .......................
        %
        % domain_structure.centers:
        % 0. it is a matrix 2 x 2;
        % 1. first row consists of the coordinates of the first center of
        %    the 'asymmetric-circular-sector';
        % 2. second row consists of the coordinates of the second center of
        %    the 'asymmetric-circular-sector';
        %
        % domain_structure.radius:
        % 0. it is a vector of length 1;
        % 1. it consists of the radius of the disk;
        %
        % domain_structure.angles:
        % 0. it is a vector of length 2;
        % 1. the first component is the parameter "alpha" and must be in
        %    [-pi,pi];
        % 2. the second component is the parameter "beta" and must be in
        %    [-pi,pi];

        % default example:
        domain_structure.centers=[1 1; 2 2];
        domain_structure.radius=3;
        domain_structure.angles=[pi/4; pi/2];

        while k <= L
            strL=args{k};
            switch strL
                case 'centers'
                    domain_structure.centers=args{k+1};
                case 'radius'
                    domain_structure.radius=args{k+1};
                case 'angles'
                    domain_structure.angles=args{k+1};
            end
            k=k+2;
        end

        domain_structure.NURBS=ellblend2NURBS(domain_structure);


        % ellblend structure

        R=domain_structure.radius;
        C=domain_structure.centers;
        thetaV=domain_structure.angles;

        domain_structure.A=[R 0; 0 0];
        domain_structure.B=[0 R; 0 0];
        domain_structure.C=C;
        domain_structure.alpha=thetaV(1);
        domain_structure.beta=thetaV(2);



    case 'asymmetric-annulus'

        % ........................ setting defaults .......................
        %
        % IMPORTANT: the 'asymmetric-annulus' is determined by a first disk
        % "D1" and a second disk "D2". It is required that "D2" is included
        % in "D1".
        %
        % domain_structure.centers:
        % 0. it is a matrix 2 x 2;
        % 1. first row consists of the coordinates of the first center of
        %    the 'asymmetric-circular-sector';
        % 2. second row consists of the coordinates of the second center of
        %    the 'asymmetric-circular-sector';
        %
        % domain_structure.radii:
        % 0. it is a vector of length 2;
        % 1. first value consists of the radius of the first disk;
        % 2. second value consists of the radius of the second disk;
        %


        % default example:
        domain_structure.centers=[0.5 0.5; 0.7 0.7];
        domain_structure.radii=[0.5 0.15];

        while k <= L
            strL=args{k};
            switch strL
                case 'centers'
                    domain_structure.centers=args{k+1};
                case 'radii'
                    domain_structure.radii=args{k+1};
            end
            k=k+2;
        end

        
        % ellblend structure
        
        R=domain_structure.radii;
        R1=R(1); R2=R(2);
        C=domain_structure.centers;

        domain_structure.A=[R1 0; R2 0];
        domain_structure.B=[0 R1; 0 R2];
        domain_structure.C=C;
        domain_structure.alpha=-pi;
        domain_structure.beta=pi;
        


    case 'vertical-circular-zone'

        % ........................ setting defaults .......................
        %
        % domain_structure.center:
        % 0. it is a matrix 1 x 2;
        % 1. it consists of the coordinates of center of the
        %    'vertical-circular-zone';

        % domain_structure.radius:
        % 0. it is a vector of length 1;
        % 1. it consists of the radius of the disk;

        % domain_structure.angles:
        % 0. it is a vector of length 2;
        % 1. the first component is the parameter "alpha" and must be in
        %    [0,pi];
        % 2. the second component is the parameter "beta" and must be in
        %    [0,pi]; it is required that "alpha < beta";

        domain_structure.center=[1 1];
        domain_structure.radius=1;
        domain_structure.angles=[pi/6 pi/2+pi/6];

        while k <= L
            strL=args{k};
            switch strL
                case 'center'
                    domain_structure.center=args{k+1};
                case 'radius'
                    domain_structure.radius=args{k+1};
                case 'angles'
                    domain_structure.angles=args{k+1};
            end
            k=k+2;
        end

        domain_structure.NURBS=ellblend2NURBS(domain_structure);

        % ellblend structure
        
        R=domain_structure.radius;
        C=domain_structure.center;
        thetaV=domain_structure.angles;

        domain_structure.A=[R 0; R 0];
        domain_structure.B=[0 R; 0 -R];
        domain_structure.C=[C; C];
        domain_structure.alpha=thetaV(1);
        domain_structure.beta=thetaV(2);



    case 'horizontal-circular-zone'
        % ........................ setting defaults .......................
        %
        % domain_structure.center:
        % 0. it is a matrix 1 x 2;
        % 1. it consists of the coordinates of center of the
        %    'vertical-circular-zone';
        %
        % domain_structure.radius:
        % 0. it is a vector of length 1;
        % 1. it consists of the radius of the disk;
        %
        % domain_structure.angles:
        % 0. it is a vector of length 2;
        % 1. the first component is the parameter "alpha" and must be in
        %    [-pi/2,pi/2];
        % 2. the second component is the parameter "beta" and must be in
        %    [-pi/2,pi/2]; it is required that "alpha < beta";

        domain_structure.center=[1 1];
        domain_structure.radius=1;
        domain_structure.angles=[-pi/4 pi/6];

        while k <= L
            strL=args{k};
            switch strL
                case 'center'
                    domain_structure.center=args{k+1};
                case 'radius'
                    domain_structure.radius=args{k+1};
                case 'angles'
                    domain_structure.angles=args{k+1};
            end
            k=k+2;
        end

        domain_structure.NURBS=ellblend2NURBS(domain_structure);

        % ellblend structure
        
        R=domain_structure.radius;
        C=domain_structure.center;
        thetaV=domain_structure.angles;

        domain_structure.A=[R 0; -R 0];
        domain_structure.B=[0 R; 0 R];
        domain_structure.C=[C; C];
        domain_structure.alpha=thetaV(1);
        domain_structure.beta=thetaV(2);



    case 'circular-segment'

        % ........................ setting defaults .......................
        %
        % domain_structure.center:
        % 0. it is a matrix 1 x 2;
        % 1. it consists of the coordinates of center of the
        %    'circular-segment';
        %
        % domain_structure.radius:
        % 0. it is a vector of length 1;
        % 1. it consists of the radius of the disk;
        %
        % domain_structure.angles:
        % 0. it is a vector of length 2;
        % 1. the first component is the parameter "alpha" and must be in
        %    [-pi,pi];
        % 2. the second component is the parameter "beta" and must be in
        %    [-pi,pi];

        domain_structure.center=[1 1];
        domain_structure.radius=2;
        domain_structure.angles=[-pi/4 pi/6];

        while k <= L
            strL=args{k};
            switch strL
                case 'center'
                    domain_structure.center=args{k+1};
                case 'radius'
                    domain_structure.radius=args{k+1};
                case 'angles'
                    domain_structure.angles=args{k+1};
            end
            k=k+2;
        end

        % NURBS structure

        domain_structure.NURBS=ellblend2NURBS(domain_structure);

        % ellblend structure
        
        R=domain_structure.radius;
        C=domain_structure.center;
        thetaV=domain_structure.angles;
        alpha=thetaV(1); beta=thetaV(2); gamma=(thetaV(1)+thetaV(2))/2;

        domain_structure.A=[R*cos(alpha) R*sin(alpha); R*cos(beta) R*sin(beta)];
        domain_structure.B=[-R*sin(alpha) R*cos(alpha); R*sin(beta) -R*cos(beta)];

        domain_structure.C=[C; C];
        domain_structure.alpha=0;
        domain_structure.beta=abs(beta-alpha)/2;



    case 'symmetric-lens'

        % ........................ setting defaults .......................
        %
        % domain_structure.center:
        % 0. it is a matrix 1 x 1;
        % 1. it consists of the x0 coordinates of center of the
        %    'symmetric-lens', i.e. the lens will be defined by the disks
        %    having centers [+/- domain_structure.center,0];
        %
        % domain_structure.radius:
        % 0. it is a vector of length 1;
        % 1. it consists of the radius of each disk of the lens;
        %

        % setting defaults
        domain_structure.center=1;
        domain_structure.radius=2.5;

        while k <= L
            strL=args{k};
            switch strL
                case 'center'
                    domain_structure.center=args{k+1};
                case 'radius'
                    domain_structure.radius=args{k+1};
            end
            k=k+2;
        end


        % ellblend structure
        
        R=domain_structure.radius;
        a=domain_structure.center;

        domain_structure.A=[R 0; -R 0];
        domain_structure.B=[0 R; 0 R];
        domain_structure.C=[-a 0; a 0];

        domain_structure.alpha=-acos(a/R);
        domain_structure.beta=acos(a/R);



    case 'butterfly'

        % ........................ setting defaults .......................
        %
        % domain_structure.center:
        % 0. it is a matrix 1 x 2;
        % 1. it consists of the coordinates of center of the
        %    'butterfly';
        %
        % domain_structure.radius:
        % 0. it is a vector of length 1;
        % 1. it consists of the radius of the disk of the
        %    butterfly;
        %
        % domain_structure.angles:
        % 0. it is a vector of length 2;
        % 1. the first component is the parameter "alpha" and must be in
        %    [-pi,pi];
        % 2. the second component is the parameter "beta" and must be in
        %    [-pi,pi];

        % setting defaults
        domain_structure.center=[1 2];
        domain_structure.radius=1;
        domain_structure.angles=[-pi/4 pi/6];


        while k <= L
            strL=args{k};
            switch strL
                case 'center'
                    domain_structure.center=args{k+1};
                case 'radius'
                    domain_structure.radius=args{k+1};
                case 'angles'
                    domain_structure.angles=args{k+1};
            end
            k=k+2;
        end

        domain_structure.NURBS=ellblend2NURBS(domain_structure);

        % ellblend structure
        
        R=domain_structure.radius;
        C=domain_structure.center;
        thetaV=domain_structure.angles;

        domain_structure.A=[R 0; -R 0];
        domain_structure.B=[0 R; 0 -R];
        domain_structure.C=[C; C];
        domain_structure.alpha=thetaV(1);
        domain_structure.beta=thetaV(2);



    case 'candy'

        % ........................ setting defaults .......................
        %
        % domain_structure.center:
        % 0. it is a matrix 1 x 1;
        % 1. it determines the value of the parameter "a"  defining
        %    the candy;
        %
        % domain_structure.radius:
        % 0. it is a vector of length 1;
        % 1. it consists of the radius "r" defining the candy;
        %
        % domain_structure.angle:
        % 0. it is a vector of length 1;
        % 1. it is the parameter "alpha" defining the candy and must be in
        %    [-pi,pi];
        %
        % IMPORTANT: it must be "-alpha > acos(a/r)"


        % setting defaults
        domain_structure.center=0.5;
        domain_structure.radius=1;
        domain_structure.angle=-pi/2;


        while k <= L
            strL=args{k};
            switch strL
                case 'center'
                    domain_structure.center=args{k+1};
                case 'radius'
                    domain_structure.radius=args{k+1};
                case 'angle'
                    domain_structure.angle=args{k+1};
            end
            k=k+2;
        end


        % ellblend structure
        
        r=domain_structure.radius;
        a=domain_structure.center;
        alpha=domain_structure.angle;

        domain_structure.A=[r 0; -r 0];
        domain_structure.B=[0 r;0 r];
        domain_structure.C=[-a 0; a 0];
        domain_structure.alpha=alpha;
        domain_structure.beta=-alpha;


    case 'NURBS'
        % ........................ setting defaults .......................

        %
        % domain_structure.geometry:
        % 0. it is a "Matlab-structure" defining the geometry of
        %    the NURBS;
        % 1. for examples on how to define this "structure"
        %    see the gallery "gallery_NURBS".

        % setting default from gallery of domains
        
        example_type=42;
        domain_structure=gallery_domains_2D(domain_str,example_type);

        while k <= L
            strL=args{k};
            switch strL
                case 'geometry'
                    domain_structure.NURBS=args{k+1};
            end
            k=k+2;
        end



    case 'union-disks'

        % ........................ setting defaults .......................
        %
        % domain_structure.centers:
        % 0. it is a matrix M x 2;
        % 1. the k-th row consist of the coordinates of the center
        %    of the k-th disk defining the domain;
        %
        % domain_structure.radii:
        % 0. it is a vector of length M;
        % 1. the k-th component consist of the radius of the k-th
        %     disk defining the domain;


        % setting defaults
        example_type=1;
        domain_structure=gallery_domains_2D(domain_str,example_type);

        while k <= L
            strL=args{k};
            switch strL
                case 'centers'
                    domain_structure.centers=args{k+1};
                case 'radii'
                    domain_structure.radii=args{k+1};
            end
            k=k+2;
        end



    case 'square'

        % reference square: [-1,1] x [-1,1]

        vertices=[-1 -1; 1 -1; 1 1; -1 1; -1 -1];
        domain_structure.vertices=vertices;

        XV=domain_structure.vertices(:,1);
        YV=domain_structure.vertices(:,2);

        domain_structure.polyshape=polyshape(XV,YV);



    case 'unit-square[0,1]x[0,1]'

        % reference square: [0,1] x [0,1]

        vertices=[0 0; 1 0; 1 1; 0 1; 0 0];
        domain_structure.vertices=vertices;

        XV=domain_structure.vertices(:,1);
        YV=domain_structure.vertices(:,2);

        domain_structure.polyshape=polyshape(XV,YV);



    case 'rectangle'

        % 0. defining the domain as dbox=[xLimits; yLimits], where
        %    xLimits=[xmin, xmax], yLimits=[ymin, ymax]
        %
        % 1. the routine will also assign the vertices of the domain.

        % setting defaults
        domain_structure.xLimits=[0.5 1];
        domain_structure.yLimits=[2 3];

        while k <= L
            strL=args{k};
            switch strL
                case 'xLimits'
                    domain_structure.xLimits=args{k+1};
                case 'yLimits'
                    domain_structure.yLimits=args{k+1};
            end
            k=k+2;
        end

        xmin=min(domain_structure.xLimits);
        xmax=max(domain_structure.xLimits);
        ymin=min(domain_structure.yLimits);
        ymax=max(domain_structure.yLimits);

        vertices=[xmin ymin; xmax ymin; xmax ymax; xmin ymax; ...
            xmin ymin];

        domain_structure.vertices=vertices;

        XV=domain_structure.vertices(:,1);
        YV=domain_structure.vertices(:,2);

        domain_structure.polyshape=polyshape(XV,YV);



    case 'triangle'

        % ........................ setting defaults .......................
        %
        % domain_structure.vertices:
        % 0. it is a matrix 3 x 2 (no requirement that the first and last
        %     vertex coincide;
        % 1. each row consists of the coordinates of a vertex of the domain
        %    in cartesian coordinates;

        % setting defaults
        vertices=[0 0; 1 0; 1 1; 0 0];

        while k <= L
            strL=args{k};
            switch strL
                case 'vertices'
                    vertices=args{k+1};
            end
            k=k+2;
        end

        if norm(vertices(1,:)-vertices(end,:)) > 0
            vertices(end+1,:)=vertices(1,:);
        end

        domain_structure.vertices=vertices;

        XV=domain_structure.vertices(:,1);
        YV=domain_structure.vertices(:,2);

        domain_structure.polyshape=polyshape(XV,YV);



    case 'polygcirc'

        % ........................ setting defaults .......................
        %
        % domain_structure.extrema:
        % 0. it is a matrix 2 x 2;
        % 1. each row consists of the coordinates of a vertex of the domain
        %    in cartesian coordinates;
        % 2. the first row is the first extremum (point from which it
        %    starts the arc of the disk);
        % 3. the second row is the second extremum (point from which it
        %    end the arc of the disk);
        %
        % domain_structure.vertices:
        % 0: it is a matrix N x 2 of polygon vertices (cart. coordinates);
        % 1: arc extrema are added to "domain_structure.vertices", by
        %    default in counterclockwise order;
        %
        % domain_structure.center:
        % 0. it is a matrix 1 x 2;
        % 1. it determines the center of the disk defining "polygcirc"
        %    domain;
        %
        % domain_structure.radius:
        % 0. it is a vector of length 1;
        % 1. it consists of the radius "r" of the disk defining "polygcirc"
        %    domain;
        %
        % domain_structure.convex:
        % 0: it is a parameter that if "0" will define a concave arc,
        %     otherwise a convex arc.

        % default example
        example_type=2;
        domain_structure=gallery_domains_2D(domain_str,example_type);

        while k <= L
            strL=args{k};
            switch strL
                case 'extrema'
                    domain_structure.extrema=args{k+1};
                case 'vertices'
                    domain_structure.vertices=args{k+1};
                case 'center'
                    domain_structure.center=args{k+1};
                case 'radius'
                    domain_structure.radius=args{k+1};
                case 'convex'
                    domain_structure.convex=args{k+1};
            end
            k=k+2;
        end


    case 'unit-simplex'

        domain_structure.vertices=[0 0; 1 0; 0 1; 0 0];

        XV=domain_structure.vertices(:,1);
        YV=domain_structure.vertices(:,2);

        domain_structure.polyshape=polyshape(XV,YV);

    case 'sphere'

        domain_structure=gallery_domains_2D(domain_str);
        %XW=cub_sphere(deg);

     case 'spherical-polygon'

        example_type=2;
        domain_structure=gallery_domains_S2(domain_str,example_type);

        sphpgon_vertices=domain_structure.vertices;
        %XW=cub_sphpgon(deg,sphpgon_vertices);


        

end

domain_structure.dbox=domain_boundingbox(domain_structure);



