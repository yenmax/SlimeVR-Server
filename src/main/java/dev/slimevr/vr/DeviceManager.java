package dev.slimevr.vr;

import dev.slimevr.VRServer;
import io.eiren.util.collections.FastList;
import org.kohsuke.github.GHCommit;
import org.kohsuke.github.GHRelease;
import org.kohsuke.github.GHRepository;
import org.kohsuke.github.GitHub;

import java.io.IOException;


public class DeviceManager {

	private final FastList<Device> devices = new FastList<>();

	public FastList<Device> getDevices() {
		return devices;
	}

	private final VRServer server;

	public DeviceManager(VRServer server) {
		this.server = server;
	}

	public Device createDevice(String name, String version, String manufacturer) {
		Device device = new Device();

		device.setCustomName(name);
		device.setFirmwareVersion(version);
		device.setManufacturer(manufacturer);

		return device;
	}

	public void addDevice(Device device) {
		this.server.queueTask(() -> {
			this.devices.add(device);
			try {
				this.canUpdateDevice(device);
			} catch (IOException e) {
				e.printStackTrace();
			}
		});
	}

	public void updateDeviceFirmware(Device device) {

	}

	public boolean canUpdateDevice(Device device) throws IOException {

		GitHub gitHub = GitHub.connectAnonymously();
		GHRepository repository = gitHub
			.getUser("SlimeVR")
			.getRepository("SlimeVR-Tracker-ESP");
		GHRelease release = repository
			.getLatestRelease();
		GHCommit commit = repository.getCommit(release.getTagName());
//		repository
//			.listArtifacts()
//			.withPageSize(1)
//			.forEach((a) -> System.out.println(a.getName()));
		return device.getBoardType() == BoardType.BOARD_SLIMEVR;
	}
}
