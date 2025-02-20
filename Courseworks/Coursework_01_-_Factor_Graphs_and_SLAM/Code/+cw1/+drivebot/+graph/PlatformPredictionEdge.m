classdef PlatformPredictionEdge < g2o.core.BaseBinaryEdge
    % PlatformPredictionEdge summary of PlatformPredictionEdge
    %
    % This class stores the factor representing the process model which
    % transforms the state from timestep k to k+1
    %
    % The process model is as follows.
    %
    % Define the rotation vector
    %
    %   M = dT * [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0;0 0 1];
    %
    % The new state is predicted from 
    %
    %   x_(k+1) = x_(k) + M * [vx;vy;theta]
    %
    % Note in this case the measurement is actually the mean of the process
    % noise. It has a value of 0. The error vector is given by
    %
    % e(x,z) = inv(M) * (x_(k+1) - x_(k))
    %
    % Note this requires estimates from two vertices - x_(k) and x_(k+1).
    % Therefore, this inherits from a binary edge. We use the convention
    % that vertex slot 1 contains x_(k) and slot 2 contains x_(k+1).
    
    properties(Access = protected)
        % The length of the time step
        dT;
    end
    
    methods(Access = public)
        function obj = PlatformPredictionEdge(dT)
            % PlatformPredictionEdge for PlatformPredictionEdge
            %
            % Syntax:
            %   obj = PlatformPredictionEdge(dT);
            %
            % Description:
            %   Creates an instance of the PlatformPredictionEdge object.
            %   This predicts the state from one timestep to the next. The
            %   length of the prediction interval is dT.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a PlatformPredictionEdge

            assert(dT >= 0);
            obj = obj@g2o.core.BaseBinaryEdge(3);            
            obj.dT = dT;
        end
       
        function initialEstimate(obj)
            % INITIALESTIMATE Compute the initial estimate of a platform.
            %
            % Syntax:
            %   obj.initialEstimate();
            %
            % Description:
            %   Compute the initial estimate of the platform x_(k+1) given
            %   an estimate of the platform at time x_(k) and the control
            %   input u_(k+1)

            % Retrieve the prior state estimate x_k.
            priorX = obj.edgeVertices{1}.estimate();
            
            % Extract the heading (yaw) angle from x_k.
            psi = priorX(3);
            
            % Compute the rotation matrix M(ψ_k).
            M = [cos(psi) -sin(psi) 0;
                 sin(psi)  cos(psi) 0;
                 0         0        1];
            
            % Compute the predicted state assuming zero process noise:
            % Multiply the control input by dT to convert from a velocity to a displacement.
            predictedX = priorX + M * (obj.dT * obj.z);
            
            % Normalize the heading angle to the range [-pi, pi].
            predictedX(3) = g2o.stuff.normalize_theta(predictedX(3));
            
            % Set the estimated state for the posterior vertex.
            obj.edgeVertices{2}.setEstimate(predictedX);
        end
        
        function computeError(obj)
            % COMPUTEERROR Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the measurement and the parameter state in the
            %   vertex. Note the error enters in a nonlinear manner, so the
            %   equation has to be rearranged to make the error the subject
            %   of the formulat
                       
            % Get the state x_k from the prior vertex.
            priorX = obj.edgeVertices{1}.x;
            psi = priorX(3);
            
            % Compute the inverse rotation matrix M(ψ_k)^{-1}.
            % Since M(ψ) is an orthonormal rotation matrix, its inverse is its transpose.
            % We write it out explicitly:
            M_inv = [cos(psi) sin(psi) 0;
                     -sin(psi) cos(psi) 0;
                     0         0        1];
            
            % Compute the state difference (x_(k+1) - x_k).
            deltaX = obj.edgeVertices{2}.x - priorX;
            
            % Compute the error: transform the state difference into the vehicle's frame
            % and subtract the expected displacement (dT scaled control input).
            error = M_inv * deltaX - obj.dT * obj.z;
            
            % Normalize the heading error.
            error(3) = g2o.stuff.normalize_theta(error(3));
            
            % Store the computed error.
            obj.errorZ = error;
        end
        
        % Compute the Jacobians
        function linearizeOplus(obj)
            % LINEARIZEOPLUS Compute the Jacobians for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the Jacobians for the edge. Since we have two
            %   vertices which contribute to the edge, the Jacobians with
            %   respect to both of them must be computed.
            %

            % Retrieve the prior state x_k.
            priorX = obj.edgeVertices{1}.x;
            psi = priorX(3);
            c = cos(psi);
            s = sin(psi);
            
            % The inverse rotation matrix M(ψ_k)^{-1}:
            M_inv = [c  s  0;
                    -s  c  0;
                     0  0  1];
            
            % The state difference (x_(k+1) - x_k).
            deltaX = obj.edgeVertices{2}.x - priorX;
            
            % Jacobian with respect to x_(k+1) is simply M(ψ_k)^{-1}.
            obj.J{2} = M_inv;
            
            % Jacobian with respect to x_k:
            % We have two contributions:
            % 1. The derivative of - (x_(k+1)-x_k) gives -I transformed by M_inv: -M_inv.
            % 2. The derivative of M_inv with respect to ψ (which only affects the third column)
            %    multiplied by (x_(k+1)-x_k). Using the derivation in the VehicleKinematicsEdge:
            %
            obj.J{1} = zeros(3,3);
            obj.J{1}(1,1) = -c;
            obj.J{1}(1,2) = -s;
            obj.J{1}(1,3) = -deltaX(1)*s + deltaX(2)*c;
            
            obj.J{1}(2,1) = s;
            obj.J{1}(2,2) = -c;
            obj.J{1}(2,3) = -deltaX(1)*c - deltaX(2)*s;
            
            % For the yaw component, the derivative is simply -1.
            obj.J{1}(3,3) = -1;
            
            % Note: The dT * u term does not contribute to the Jacobians because it is constant.
        end
    end    
end