function [PC, RC] = ssa(data, showVals)
    M = 70;    % window length = embedding dimension
    N = length(data);   % length of generated time series
    t = (1:N)';
    % normalise the data
    data = data - mean(data);            % remove mean value
    data = data/std(data,1);             % normalize to standard deviation 1

    % cov matric C (trajectory)
    % can't remember if this is the best method to use
    % believe Jarchi paper uses it but not other methods
    Y=zeros(N-M+1,M);
    for m=1:M
      Y(:,m) = data((1:N-M+1)+m-1);
    end
    Cemb=Y'*Y / (N-M+1);
    C = Cemb;

    % calculate LAMBDA (eigenvals) and RHO (eigenvectors)
    [RHO,LAMBDA] = eig(C);
    LAMBDA = diag(LAMBDA);               % extract the diagonal elements
    [LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
    RHO = RHO(:,ind);                    % and eigenvectors

    % show the eigenvals and eigenvectors (optional)
    if showVals
        figure;
        set(gcf,'name','Eigenvectors RHO and eigenvalues LAMBDA')
        clf;
        subplot(3,1,1);
        plot(LAMBDA,'o-');
        subplot(3,1,2);
        plot(RHO(:,1:2), '-');
        legend('1', '2');
        subplot(3,1,3);
        plot(RHO(:,3:4), '-');
        legend('3', '4');
    end

    % calculate PCs
    PC = Y*RHO;

    % show the PCs
    if showVals
        figure;
        set(gcf,'name','Principal components PCs')
        clf;
        for m=1:4
          subplot(4,1,m);
        %   plot(t(1:N-M+1), PC(:,m),'k-');
          plot(PC(:,m),'k-');
          ylabel(sprintf('PC %d',m));
          ylim([-10 10]);
        end
    end

    % calculate RCs
    RC=zeros(N,M);
    for m=1:M
      buf=PC(:,m)*RHO(:,m)'; % invert projection
      buf=buf(end:-1:1,:);
      for n=1:N % anti-diagonal averaging
        RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
      end
    end

    % show the RCs
    if showVals
        figure;
        set(gcf,'name','Reconstructed components RCs')
        clf;
        for m=1:4
          subplot(4,1,m);
          plot(t,RC(:,m),'r-');
          ylabel(sprintf('RC %d',m));
          ylim([-1 1]);
        end
    end

    % show the reconstructed signals
    figure;
    set(gcf,'name','Original time series data and reconstruction RC')
    title('Original ML Plane Signals vs Reconstruction')
    clf;
    subplot(4,1,1)
    plot(t,data,'b-',t,sum(RC(:,:),2),'r-');
    legend('Original','All RCs');
    
    subplot(4,1,2)
    plot(t,data,'b','LineWidth',2);
    plot(t,data,'b-',t,sum(RC(:,1:2),2),'r-');
    legend('Original','RCs 1-2');
    
    subplot(4,1,3)
    plot(t,data,'b','LineWidth',2);
    plot(t,data,'b-',t,sum(RC(:,1:4),2),'r-');
    legend('Original','1-4');
    
    subplot(4,1,4)
    plot(t,data,'b','LineWidth',2);
    plot(t,data,'b-',t,sum(RC(:,5:end),2),'r-');
    legend('Original','5-end');

end





