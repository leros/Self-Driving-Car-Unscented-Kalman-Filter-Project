[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

---
### PART 2: Sensor Fusion, Localization, and Control
### Project 2:  Unscented Kalman Filter Project

I implement an unscented Kalman filter(UKF) using the CTRV motion model. I use the same bicycle simulation data set from the extended Kalman filter(EKF) project. That way I can compare your results with the EKF project.

A standard Kalman filter can only handle linear equations. Both the EKF and the UKF allow people to use non-linear equations; the difference between EKF and UKF is how they handle non-linear equations. But the basics are the same: initialize, predict, update.


Passing the project requires obtaining RMSE values that are lower that the tolerance outlined in the [project rubric](https://review.udacity.com/#!/rubrics/783/view). 