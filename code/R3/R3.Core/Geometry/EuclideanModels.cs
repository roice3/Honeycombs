namespace R3.Geometry
{
	using R3.Geometry;

	public enum EuclideanModel
	{
		Isometric,
		Conformal,
		Disk,
		UpperHalfPlane,
	}

	public class EuclideanModels
	{
		public static Vector3D DiskToIsometric( Vector3D v )
		{
			// ZZZ - Check that this is correct (it's quite close if not!)
			return SphericalModels.StereoToGnomonic( v );
		}

		public static Vector3D UpperHalfPlaneToIsometric( Vector3D v )
		{
			v = HyperbolicModels.UpperToPoincare( v );
			v = SphericalModels.StereoToGnomonic( v );
			return v;	
		}
	}
}
