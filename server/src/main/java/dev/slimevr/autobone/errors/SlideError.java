package dev.slimevr.autobone.errors;

import dev.slimevr.autobone.AutoBoneTrainingStep;
import dev.slimevr.tracking.processor.skeleton.HumanSkeleton;
import dev.slimevr.tracking.trackers.ComputedTracker;
import dev.slimevr.tracking.trackers.TrackerRole;


// The change in position of the ankle over time
public class SlideError implements IAutoBoneError {
	@Override
	public float getStepError(AutoBoneTrainingStep trainingStep) throws AutoBoneException {
		return getSlideError(
			trainingStep.getHumanPoseManager1().getSkeleton(),
			trainingStep.getHumanPoseManager2().getSkeleton()
		);
	}

	public static float getSlideError(HumanSkeleton skeleton1, HumanSkeleton skeleton2) {
		// Calculate and average between both feet
		return (getSlideError(skeleton1, skeleton2, TrackerRole.LEFT_FOOT)
			+ getSlideError(skeleton1, skeleton2, TrackerRole.RIGHT_FOOT)) / 2f;
	}

	public static float getSlideError(
		HumanSkeleton skeleton1,
		HumanSkeleton skeleton2,
		TrackerRole trackerRole
	) {
		// Calculate and average between both feet
		return getSlideError(
			skeleton1.getComputedTracker(trackerRole),
			skeleton2.getComputedTracker(trackerRole)
		);
	}

	public static float getSlideError(ComputedTracker tracker1, ComputedTracker tracker2) {
		// Return the midpoint distance
		return tracker1.position.distance(tracker2.position) / 2f;
	}
}
