clear all;
N = 81;
Radius=1;
NPixels = 9;
NBLocks= N/NPixels;
tmd = round(rand(N*N));
tm2d = round(rand(N*N));

    for i = 1:NPixels:N+1
        for j= 1:NPixels:N+1
            for dx=-Radius:1:Radius
                for dy=-Radius:1:Radius
                    cond=(i + dx*NPixels >= 0) & (i + dx*NPixels < N);
%                     and (j + dy*NPixels >= 0) and (j + dy*NPixels < N);
%                  if ()
%                      
%                     temp = temp + abs(tmd[(i+k)*N+(j+l)]-tm2d[(i+k+dx*Radius)*N+(j+l+dy*Radius)]);
%                  end
             end
        end
    end
    end