import { RenameItem } from './Item'

export const RenameList = props => (
  <div>
    <div>
      <button type="button" onClick={props.handleAdd}>Add</button>
    </div>
    <ul>
      {props.renames.map((rename, index) => (
        <RenameItem
          key={rename.uuid}
          data={rename}
          rule={props.rule}
          domains={props.domains}
          handleChange={props.handleChange(index)}
          handleRemove={props.handleRemove(index)}
        />
      ))}
    </ul>
  </div>
)
