package dev.slimevr.vr;

public enum BoardType {

	UNKNOWN(0),
	BOARD_SLIMEVR_LEGACY(1),
	BOARD_SLIMEVR_DEV(2),
	BOARD_NODEMCU(3),
	BOARD_CUSTOM(4),
	BOARD_WROOM32(5),
	BOARD_WEMOSD1MINI(6),
	BOARD_TTGO_TBASE(7),
	BOARD_ESP01(8),
	BOARD_SLIMEVR(9),
	BOARD_LOLIN_C3_MINI(10);

	private final int id;

	BoardType(int id) {
		this.id = id;
	}

	public static BoardType getBoardType(int id) {
		if (id >= BoardType.values().length || id <= 0)
			return BoardType.UNKNOWN;
		else
			return BoardType.values()[id];
	}
}
