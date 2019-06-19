bl_info = {
    "name": "Convex Hull Generation in Julia",
    "author": "Jackson Meade (NCSSM) and Jack Snoeyink (UNC Chapel Hill)",
    "version": (1, 0),
    "blender": (2, 80, 0),
    "location": "View3D > Add > Mesh > CGinJulia Convex Hull",
    "description": "Adds a new Mesh Object",
    "warning": "Please make sure you are familiar with the CGinJulia package before using this EXPERIMENTAL add-on",
    "wiki_url": "github.com/snoeyink/CGinJulia",
    "category": "Add Convex Hull",
}


import bpy
from bpy.types import Operator
from bpy.props import FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector


def add_object(self, context):
    scale_x = self.scale.x
    scale_y = self.scale.y
    scale_z = self.scale.z

    verts = [Vector((10.307114 * scale_x,59.001675 * scale_y,1.0076237 * scale_z)),Vector((50.291855 * scale_x,54.045673 * scale_y,18.409445 * scale_z)),Vector((28.580631 * scale_x,1.3269854 * scale_y,22.894176 * scale_z)),Vector((59.29259 * scale_x,30.974808 * scale_y,21.623482 * scale_z)),Vector((54.506737 * scale_x,36.03917 * scale_y,59.34175 * scale_z)),Vector((29.455175 * scale_x,10.385149 * scale_y,0.40928364 * scale_z)),Vector((26.764748 * scale_x,56.28667 * scale_y,40.8207 * scale_z)),Vector((53.86461 * scale_x,27.538862 * scale_y,24.953556 * scale_z)),Vector((28.633432 * scale_x,59.284702 * scale_y,31.155682 * scale_z)),Vector((20.903978 * scale_x,28.97039 * scale_y,54.883198 * scale_z)),Vector((4.821453 * scale_x,49.737373 * scale_y,24.635317 * scale_z)),Vector((54.06224 * scale_x,57.36933 * scale_y,56.077423 * scale_z)),Vector((23.656296 * scale_x,37.11685 * scale_y,54.310833 * scale_z)),Vector((50.47258 * scale_x,26.506662 * scale_y,15.864965 * scale_z)),Vector((14.123511 * scale_x,6.405494 * scale_y,37.599205 * scale_z)),Vector((21.397205 * scale_x,58.206276 * scale_y,44.263393 * scale_z)),Vector((41.80298 * scale_x,8.873413 * scale_y,20.180246 * scale_z)),Vector((38.05564 * scale_x,50.938583 * scale_y,53.06063 * scale_z)),]
    edges = []
    faces = [[1,6,2],[1,1,1],[2,1,6],[2,4,1],[2,1,9],[3,6,1],[3,1,5],[3,1,6],[4,6,1],[4,1,5],[5,4,1],[5,1,1],[5,1,3],[6,1,1],[6,2,1],[6,3,1],[6,4,2],[6,1,4],[9,2,1],[1,5,1],[1,1,5],[1,1,1],[1,1,1],[1,1,1],[1,1,1],[1,2,4],[1,1,1],[1,5,1],[1,1,1],[1,1,1],[1,6,1],[1,1,1],[1,9,1],[1,1,1],[1,1,1],[1,5,4],[1,6,3],]
    mesh = bpy.data.meshes.new(name="Convex Hull of HelloWorld.py")
    mesh.from_pydata(verts, edges, faces)
    # useful for development when the mesh may be invalid.
    # mesh.validate(verbose=True)
    object_data_add(context, mesh, operator=self)

class OBJECT_OT_add_object(Operator, AddObjectHelper):
    """Generate the Mesh Object from Points and CGinJulia-Constructed Faces"""
    bl_idname = "mesh.add_object"
    bl_label = "Show Representation"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):

        add_object(self, context)

        return {'FINISHED'}


# Registration

def add_object_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_object.bl_idname,
        text="CGinJulia Convex Hull (HelloWorld.py)",
        icon='SURFACE_NSPHERE')

# This allows you to right click on a button and link to the manual
def add_object_manual_map():
    url_manual_prefix = "https://docs.blender.org/manual/en/dev/"
    url_manual_mapping = (
        ("bpy.ops.mesh.add_object", "editors/3dview/object"),
    )
    return url_manual_prefix, url_manual_mapping


def register():
    bpy.utils.register_class(OBJECT_OT_add_object)
    bpy.utils.register_manual_map(add_object_manual_map)
    bpy.types.VIEW3D_MT_mesh_add.append(add_object_button)


def unregister():
    bpy.utils.unregister_class(OBJECT_OT_add_object)
    bpy.utils.unregister_manual_map(add_object_manual_map)
    bpy.types.VIEW3D_MT_mesh_add.remove(add_object_button)


if __name__ == "__main__":
    register()
