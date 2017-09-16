namespace R3.Geometry
{
	using Math = System.Math;
	using System.Diagnostics;
	using System.Linq;
	using R3.Core;

	public class VectorND
	{
		public VectorND( int dimension )
		{
			X = new double[dimension];
		}

		public VectorND( double[] components )
		{
			X = components;
		}

		public Vector3D ToVec3D()
		{
			return new Vector3D( X[0], X[1], X[2], X[3] );
		}

		public VectorND Clone()
		{
			return new VectorND( (double[])this.X.Clone() );
		}

		public int Dimension 
		{
			get { return X.Length; }
			set { X = new double[value]; }
		}

		public double[] X { get; set; }

		public static VectorND operator /( VectorND v, double s )
		{
			double[] components = new double[v.Dimension];

			for( int i = 0; i < components.Length; i++ )
				components[i] = v.X[i] / s;

			return new VectorND( components );
		}

		public void Divide( double s )
		{
			for( int i = 0; i < Dimension; i++ )
				X[i] /= s;
		}

		public static VectorND operator *( VectorND v, double s )
		{
			double[] components = new double[v.Dimension];

			for( int i = 0; i < components.Length; i++ )
				components[i] = v.X[i] * s;

			return new VectorND( components );
		}

		public static VectorND operator *( double s, VectorND v )
		{
			return v * s;
		}

		public static VectorND operator +( VectorND v1, VectorND v2 )
		{
			Debug.Assert( v1.Dimension == v2.Dimension );
			double[] components = new double[v1.Dimension];

			for( int i = 0; i < components.Length; i++ )
				components[i] = v1.X[i] + v2.X[i];

			return new VectorND( components );
		}

		public static VectorND operator -( VectorND v )
		{
			double[] components = new double[v.Dimension];

			for( int i = 0; i < components.Length; i++ )
				components[i] = -v.X[i];

			return new VectorND( components );
		}

		public static VectorND operator -( VectorND v1, VectorND v2 )
		{
			return v1 + ( -v2 );
		}

		public double Dot( VectorND v )
		{
			double dot = 0;
			for( int i = 0; i < this.Dimension; i++ )
				dot += this.X[i] * v.X[i];
			return dot;
		}
	}
}
