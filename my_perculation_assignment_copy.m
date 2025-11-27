fprintf('start calculate the density of spanning cluasters...\n\n');

L = 10;
p = 0.6;
fprintf('size of grid: %dx%d\n\n',L,L);
fprintf('probability of occupied grid" %.2f\n\n',p);

fprintf('1.build random grid...\n');
grid = rand(L,L)<p;

fprintf('2.present grid : \n');
disp('represent 1 = (occupate, 0 = empty):');
disp(grid);

figure(1);
imagesc(grid);
colormap([1,1,1;0,0,0]);
title('random grid(white is empty,black is occupied)');
colorbar;
axis equal;

fprintf('3.check the connection...\n');

visited = false(L,L);
found_spanning = false;
spanning_cluster_size = 0; 

for j = 1:L
    if grid(1,j)== 1&&~visited(1,j)
        fprintf('start to search form (1,%d)',j);

        stack = [1,j];
        current_cluster_size = 0;
        reaches_bottom = false;
        temp_visited = visited;

        while ~isempty(stack)
            current = stack(end,:);
            stack = stack(1:end-1,:);
            x = current(1);y=current(2);

            if ~temp_visited(x,y)

                temp_visited(x,y) = true;
                current_cluster_size = current_cluster_size +1;

                if x == L
                    reaches_bottom = true;
                end

            

               if x>1 && grid(x-1,y)==1&&~temp_visited(x-1,y)
                   stack = [stack;x-1,y];
               end

               if x<L && grid(x+1,y)==1&&~temp_visited(x+1,y)
                   stack = [stack;x+1,y];
               end

               if y>1 && grid(x,y-1)==1&&~temp_visited(x,y-1)
                   stack = [stack;x,y-1];
               end

               if y<L && grid(x,y+1)==1&&~temp_visited(x,y+1)
                   stack = [stack;x,y+1];
               end

            end
        end


        if reaches_bottom
            found_spanning = true;
            spanning_cluster_size = current_cluster_size;
            visited = temp_visited;
            fprintf('found spanning clusters! size: %d site\n',spanning_cluster_size);
            break;
        end
    end
end



fprintf('\n4.calculation result:\n');
if found_spanning
    total_sites = L*L;
    density = spanning_cluster_size/total_sites;
    fprintf('spanning cluster did exist!\n');
    fprintf('the size of it is : %d site',spanning_cluster_size);
    fprintf('the total sites of grid are : %d sites',total_sites);
    fprintf('the density of the spanning clusters are : %d',density);

    
    figure(2);
    imagesc(visited);
    colormap([1,1,1;1,0,0]);
    title('spanning cluster present(red = spanning cluster site)');
    colorbar;
    axis equal;
else
    fprintf('spanning cluster unfounded\n');
    fprintf('the density of the spanning cluster P(%.2f,%d) = 0\n',p,L);
end

fprintf('\n all donn!\n')


            
                
             

