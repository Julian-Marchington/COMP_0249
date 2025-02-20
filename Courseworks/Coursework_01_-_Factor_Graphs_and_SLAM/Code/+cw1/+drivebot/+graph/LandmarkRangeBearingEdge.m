classdef LandmarkRangeBearingEdge < g2o.core.BaseBinaryEdge
    % LandmarkRangeBearingEdge summary of LandmarkRangeBearingEdge
    %
    % This class stores an edge which represents the factor for observing
    % the range and bearing of a landmark from the vehicle. Note that the
    % sensor is fixed to the platform.
    %
    % The measurement model is
    %
    %    z_(k+1)=h[x_(k+1)]+w_(k+1)
    %
    % The measurements are r_(k+1) and beta_(k+1) and are given as follows.
    % The sensor is at (lx, ly).
    %
    %    dx = lx - x_(k+1); dy = ly - y_(k+1)
    %
    %    r(k+1) = sqrt(dx^2+dy^2)
    %    beta(k+1) = atan2(dy, dx) - theta_(k+1)
    %
    % The error term
    %    e(x,z) = z(k+1) - h[x(k+1)]
    %
    % However, remember that angle wrapping is required, so you will need
    % to handle this appropriately in compute error.
    %
    % Note this requires estimates from two vertices - x_(k+1) and l_(k+1).
    % Therefore, this inherits from a binary edge. We use the convention
    % that vertex slot 1 contains x_(k+1) and slot 2 contains l_(k+1).
    
    methods(Access = public)
    
        function obj = LandmarkRangeBearingEdge()
            % LandmarkRangeBearingEdge for LandmarkRangeBearingEdge
            %
            % Syntax:
            %   obj = LandmarkRangeBearingEdge(landmark);
            %
            % Description:
            %   Creates an instance of the LandmarkRangeBearingEdge object.
            %   Note we feed in to the constructor the landmark position.
            %   This is to show there is another way to implement this
            %   functionality from the range bearing edge from activity 3.
            %
            % Inputs:
            %   landmark - (2x1 double vector)
            %       The (lx,ly) position of the landmark
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a ObjectGPSMeasurementEdge

            obj = obj@g2o.core.BaseBinaryEdge(2);
        end
        
        function initialEstimate(obj)
            % INITIALESTIMATE Compute the initial estimate of the landmark.
            %
            % Syntax:
            %   obj.initialEstimate();
            %
            % Description:
            %   Compute the initial estimate of the landmark given the
            %   platform pose and observation.

            % Retrieve the current vehicle pose (state x) from vertex slot 1.
            vehicle = obj.edgeVertices{1}.x;  % Expected to be a 3x1 vector: [x; y; theta]
            
            % Extract the range and bearing measurements.
            r = obj.z(1);
            beta = obj.z(2);
            
            % Compute the landmark's position in global coordinates.
            theta = vehicle(3);
            landmark_estimate = vehicle(1:2) + r * [cos(theta + beta); sin(theta + beta)];
            
            % Set the estimate for the landmark vertex.
            obj.edgeVertices{2}.setEstimate(landmark_estimate);
        end
        
        function computeError(obj)
            % COMPUTEERROR Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the predicted and actual range-bearing measurement.

            % Retrieve the vehicle pose and landmark position.
            vehicle = obj.edgeVertices{1}.x;  % [x; y; theta]
            landmark = obj.edgeVertices{2}.x; % [l_x; l_y]
            
            % Compute the differences.
            dx = landmark(1) - vehicle(1);
            dy = landmark(2) - vehicle(2);
            
            % Predicted range.
            r_pred = sqrt(dx^2 + dy^2);
            
            % Predicted bearing.
            beta_pred = atan2(dy, dx) - vehicle(3);
            beta_pred = g2o.stuff.normalize_theta(beta_pred);
            
            % Assemble the predicted measurement.
            h = [r_pred; beta_pred];
            
            % Compute the error: measurement - prediction.
            error = obj.z - h;
            % Normalize the bearing component of the error.
            error(2) = g2o.stuff.normalize_theta(error(2));
            
            % Store the error.
            obj.errorZ = error;
        end
        
        function linearizeOplus(obj)
            % linearizeOplus Compute the Jacobian of the error in the edge.
            %
            % Syntax:
            %   obj.linearizeOplus();
            %
            % Description:
            %   Compute the Jacobian of the error function with respect to
            %   the vertex.
            %

            % Retrieve the vehicle state and landmark estimate.
            vehicle = obj.edgeVertices{1}.x;  % [x; y; theta]
            landmark = obj.edgeVertices{2}.x; % [l_x; l_y]
            
            % Compute differences.
            dx = landmark(1) - vehicle(1);
            dy = landmark(2) - vehicle(2);
            r = sqrt(dx^2 + dy^2);
            
            % To avoid division by zero, one might check if r is very small.
            % Here we assume r > 0.
            
            % Jacobian with respect to the vehicle pose (2x3 matrix).
            J_vehicle = zeros(2, 3);
            J_vehicle(1,1) = dx / r;
            J_vehicle(1,2) = dy / r;
            J_vehicle(1,3) = 0;
            
            J_vehicle(2,1) = -dy / (r^2);
            J_vehicle(2,2) = dx / (r^2);
            J_vehicle(2,3) = 1;
            
            obj.J{1} = J_vehicle;
            
            % Jacobian with respect to the landmark position (2x2 matrix).
            J_landmark = zeros(2, 2);
            J_landmark(1,1) = -dx / r;
            J_landmark(1,2) = -dy / r;
            
            J_landmark(2,1) = dy / (r^2);
            J_landmark(2,2) = -dx / (r^2);
            
            obj.J{2} = J_landmark;
        end  
    end
end