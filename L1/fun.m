function [FDx,FDy,BDx,BDy] = fun

FDx = @(U) ForwardDx(U);
FDy = @(U) ForwardDy(U);
BDx = @(U)BackwardDx(U);
BDy = @(U)BackwardDy(U);


    function Dux = ForwardDx(U)
    Dux = [U(:,2:end)-U(:,1:end-1), U(:,1) - U(:,end)];
    end

    function Duy = ForwardDy(U)
    Duy = [U(2:end,:)-U(1:end-1,:); U(1,:) - U(end,:)];
    end

    function BDux = BackwardDx(U)
        BDux = [U(:,1) - U(:,end), U(:,2:end)-U(:,1:end-1)];
    end

    function BDuy = BackwardDy(U)
        BDuy = [U(1,:) - U(end,:); U(2:end,:)-U(1:end-1,:)];
    end
        
end

