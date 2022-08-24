function [ ] = display_clust(X, G, cmap)
%%

switch_out = 1;
smb = 'o';
col = [];
% if a segmentation is provided
num_clust = max(G);
hold all;

for i =1:max(G)
    id = G==i;
    x = X(1,id);
    y = X(2,id);
    if(isempty(col))
        % scatter(x,y,80,cmap(i,:),'filled','Marker',smb,'MarkerEdgeColor',cmap(i,:));
        scatter(x,y,10,cmap(i,:),'filled','Marker',smb,'MarkerEdgeColor','k');
    else
        scatter(x,y,8,col,'Marker',smb,'MarkerEdgeColor',col);
    end
    if(switch_out==1)
        scatter(X(1,G==0),X(2,G==0),8,[0.3,0.3,0.3],'filled','Marker',smb,'MarkerEdgeColor',[0.2,0.2,0.2]);
    end
end

