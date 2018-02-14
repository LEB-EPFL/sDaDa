function str=getSurfAndVolFromCellInfo(str,smooth)
    % calculates geometrical properties of detected meshes in structure
    % "str", corresponding to a single cell in "cellList", recording to the
    % same structure. The properties are: "steplength", "steparea",
    % "stepvolume" (length, area, and volume of each segment), "length",
    % "area", "volume" (for the whole cell), "lengthvector" - coordinates
    % along cell centerline of the centers of each segment.
    % To run this for all cells and frames of cellInfo, use getSAVforCellInfo
    if isempty(str), return; end
    if isfield(str,'mesh') && length(str.mesh)>1
        mesh = str.mesh;
        lng = size(mesh,1)-1;
        if ~isfield(str,'polarity'), str.polarity=0; end

        %length
        str.meshSteplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
            mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
        str.length = sum(str.meshSteplength);

        str.lengthvector = cumsum(str.meshSteplength)-str.meshSteplength/2;

        % area of image (2D projection)
        str.steparea = [];
        for i=1:lng, str.steparea=[str.steparea;polyarea([mesh(i:i+1,1);mesh(i+1:-1:i,3)],[mesh(i:i+1,2);mesh(i+1:-1:i,4)])]; end
        str.meshArea = sum(str.steparea);
        
        d = edist(mesh(:,1),mesh(:,2),mesh(:,3),mesh(:,4));
        % Volume and surface area, assuming radial symmetry
        [str.meshStepSA, str.meshStepvolume]=compute3DShape2([0; str.lengthvector],d);
        str.meshVolume = str.meshStepvolume(end);
        str.meshSA = str.meshStepSA(end);
        
    elseif isfield(str,'contour') && length(str.contour)>1
        contour = str.contour;
        lng = size(contour,1);
        str.length = sqrt(max(max( (repmat(contour(:,1),1,lng)-repmat(contour(:,1)',lng,1)).^2 + ...
                                   (repmat(contour(:,2),1,lng)-repmat(contour(:,2)',lng,1)).^2)));
        str.meshArea = polyarea(contour(:,1),contour(:,2));
    end
    if isfield(str,'diamFWHM')
        if max(str.diamFWHM(2,:)) < 2 % filter out diamFWHM that probably came from wrong figure (FtsZfig usually)
            str=rmfield(str,'diamFWHM');
            str=rmfield(str,'fwhmSurf');
            str=rmfield(str,'fwhmVol');
        end
        try % Volume and surface area, assuming radial symmetry
            if smooth
                pp=splinefit(str.diamFWHM(1,:),str.diamFWHM(2,:),20); % 20 breakpoints seems good.
                ySmth=ppval(pp,str.diamFWHM(1,:));
                [fwhmSurface, fwhmVolume]=compute3DShape2(str.diamFWHM(1,:)',(ySmth./2)');
            else
                [fwhmSurface, fwhmVolume]=compute3DShape2(str.diamFWHM(1,:)',(str.diamFWHM(2,:)./2)'); % /2 because we need radius not diameter!
            end
            str.fwhmSurf=fwhmSurface(end);
            str.fwhmVol=fwhmVolume(end);
        catch            
        end
    end
end

function d=edist(x1,y1,x2,y2)
    % complementary for "getextradata", computes the length between 2 points
    d=sqrt((x2-x1).^2+(y2-y1).^2);
end